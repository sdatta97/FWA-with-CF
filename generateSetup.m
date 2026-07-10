function [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,pilotIndex,D_FWA,D_cell,APpositions,UEpositions,distances,distancesUEs,isIndoor,hUT,O2IdB,repDonorGaindB] = generateSetup(params,seed)
%Generates one realization of the SMa simulation setup: sectored BS
%antenna gains (TR 38.901 Table 7.3-1), SMa pathloss / LOS probability /
%O2I penetration (TR 38.901 Clause 7.4), correlated shadow fading,
%spatial correlation matrices, repeater antenna gains, and the
%sector-association matrices D_FWA / D_cell.
%Adapted from the cell-free setup generator accompanying:
%O. T. Demir, E. Bjornson, L. Sanguinetti, "Foundations of User-Centric
%Cell-Free Massive MIMO", Foundations and Trends in Signal Processing,
%2021. Original code licensed under GPLv2.
M_sectors = size(params.locationsBS,1);
K_FWA = params.numCPE;
K = M_sectors*params.numUE+params.numCPE;
Lmax = params.Lmax;
N = params.num_antennas_per_gNB;
N_UE_FWA = params.N_UE_FWA;
N_UE_cell = params.N_UE_cell;
coverageRange_sub6 = params.coverageRange_sub6;
%% Define simulation setup

%Set the seed number if it is specified other than zero
if (nargin>1)&&(seed>0)
    rng(seed)
end

%Communication bandwidth (Hz)
B = params.Band;

%Noise figure (in dB)
noiseFigure = params.noiseFigure;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters of the log-distance model in (5.42), used only for
%the inter-UE (CPE repeater / UE-UE) links; BS-UT links use the SMa
%models of TR 38.901 Clause 7.4 below
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Sector antenna parameters (element pattern of 3GPP TR 38.901 Table 7.3-1)
S = params.sectors_per_site;             %sectors per site (co-located BS entries)
M_sites = M_sectors/S;
boresights = params.sector_boresights;   %boresight azimuth of each BS entry (deg)
tilt_zenith = params.tilt_zenith;        %mechanical downtilt, zenith convention (deg)
theta_3dB = params.theta_3dB;
phi_3dB = params.phi_3dB;
SLA_V = params.SLA_V;
A_max = params.A_max;
G_ant_max = params.G_ant_max;

%SMa propagation parameters (TR 38.901 Clause 7.4)
c0 = 3e8;
fc_GHz = params.fc/1e9;      %centre frequency (GHz); SMa PL valid for 0.5-37 GHz
hBS = params.ht_bs;
h_bldg = params.h_bldg;      %average building height h (m)
W_street = params.W_street;  %average street width W (m)
%LOS probability clutter parameters (Table 7.4.2-1, SMa row)
d_clutter = 30;
h_com = 20; h_res = 8; h_veg = 15;
r_com = 0.02; r_res = 0.18; r_veg = params.r_vegetation;
%O2I building penetration, low-loss residential model (Clause 7.4.3.1)
L_glass = 2 + 0.2*fc_GHz;
L_concrete = 5 + 4*fc_GHz;
PL_tw = 5 - 10*log10(0.3*10^(-L_glass/10) + 0.7*10^(-L_concrete/10));
sigma_P = 4.4;
%SMa LOS pathloss below the breakpoint (Table 7.4.1-1); PL2 above the
%breakpoint is PL1 evaluated at d_BP plus 40log10(d3D/d_BP)
PL1_SMa = @(d3D) 20*log10(40*pi*d3D*fc_GHz/3) + min(0.03*h_bldg^1.72,10)*log10(d3D) ...
                 - min(0.044*h_bldg^1.72,14.77) + 0.002*log10(h_bldg)*d3D;

%Prepare to save results
gainOverNoisedB = zeros(M_sectors,K);
gainOverNoisedB_ue = zeros(K,K);
R_gNB = zeros(N,N,M_sectors,K);
R_cpe = zeros(N_UE_FWA,N_UE_FWA,M_sectors,K_FWA);
%Only CPEs (indices 1:K_FWA) act as inter-UE transmitters, so the third
%dimension is K_FWA rather than K (a K x K allocation would waste GBs at
%the SMa UE counts)
R_interue = zeros(N_UE_FWA,N_UE_FWA,K_FWA,K);
R_ue = zeros(N_UE_cell,N_UE_cell,M_sectors,K-K_FWA);
distances = zeros(M_sectors,K);
distancesUEs = zeros(K,K);
pilotIndex = zeros(K,1);
D_FWA = ones(M_sectors,K_FWA);
D_cell = ones(M_sectors,K-K_FWA);
    
%Sector BS positions (co-located per site) and user positions
locationsBS = params.locationsBS;
APpositions = locationsBS(:,1) + 1i*locationsBS(:,2);
%Prepare to compute UE locations  
CPE_locations = params.CPE_locations;
UE_locations = params.UE_locations;
UEpositions = [CPE_locations; UE_locations];
UEpositions = UEpositions(:,1) + 1i*UEpositions(:,2);
%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-coverageRange_sub6 0 coverageRange_sub6],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[M_sectors 1]);

%Prepare to store shadowing correlation matrix. The BS-side field is
%normalized to unit variance and scaled per link by the SMa shadow
%fading std (4/6 dB LOS below/above breakpoint, 8 dB NLOS)
shadowCorrMatrix = ones(K,K);
shadowAPrealizations = zeros(K,M_sectors);

shadowCorrMatrix_ue = sigma_sf^2*ones(K,K);
shadowCPErealizations_ue = zeros(K,K);

%Indoor state, height, and O2I penetration loss of each user. CPEs are
%outdoor-mounted at params.hr_cpe. Cellular UEs are indoors with
%probability params.indoor_UT_ratio (Table 7.2-5), in residential
%buildings with floor heights 1.5 or 4.5 m, and receive a UT-specific
%O2I loss (low-loss residential model, Clause 7.4.3.1) applied to all
%their links
isIndoor = false(K,1);
if isfield(params,'isIndoorUE') && ~isempty(params.isIndoorUE)
    %Location-determined indoor state supplied by the caller: UEs inside
    %the housing complex are indoors, UEs outside it are outdoors
    isIndoor(K_FWA+1:K) = params.isIndoorUE;
else
    %Fall back to an i.i.d. draw with the Table 7.2-5 indoor UT ratio
    isIndoor(K_FWA+1:K) = rand(K-K_FWA,1) < params.indoor_UT_ratio;
end
hUT = params.hr*ones(K,1);
hUT(1:K_FWA) = params.hr_cpe;
O2IdB = zeros(K,1);
for k = K_FWA+1:K
    if isIndoor(k)
        if rand < 0.5
            hUT(k) = 4.5;
        else
            hUT(k) = 1.5;
        end
        d2Din = min(10*rand, 10*rand); %residential depth: min of two U(0,10) m
        O2IdB(k) = PL_tw + 0.5*d2Din + sigma_P*randn;
    end
end

%Sectored antennas at the CPEs / network-controlled repeaters: every CPE
%(each one is a potential repeater) carries S sector panels with the same
%element pattern as the gNB sectors (Table 7.3-1), boresight at the
%horizon, at a random installation azimuth. A link goes through the
%best-pointing sector. Donor-side gains (CPE towards each sector BS) are
%returned in repDonorGaindB; service-side gains (CPE towards each user)
%are embedded in gainOverNoisedB_ue below.
repOrientations = 360*rand(K_FWA,1); %installation azimuth per CPE (deg)
repDonorGaindB = zeros(M_sectors,K_FWA);

%Add UEs
for k = 1:K
    
    UEposition = UEpositions(k);
    %Compute 2D distances to all sector BSs (wrap-around aware)
    [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
    %Height difference between an AP and a UE (in meters)
    distanceVertical = params.ht_bs - hUT(k);
    distances(:,k) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
    distancesUEs(:,k) = abs(UEpositions - UEposition);
    %Sector antenna gain (dB) of each sector BS towards UE k, using the
    %element radiation pattern of TR 38.901 Table 7.3-1
    sectorGaindB = zeros(M_sectors,1);
    for l = 1:M_sectors
        %azimuth of UE k relative to sector boresight, wrapped to [-180,180]
        azimuthToUE = rad2deg(angle(UEposition - APpositionsWrapped(l,whichpos(l))));
        relAzimuth = mod(azimuthToUE - boresights(l) + 180, 360) - 180;
        %zenith angle of the LOS direction at the BS (90 deg = horizon)
        zenithToUE = 90 + rad2deg(atan(distanceVertical/max(distanceAPstoUE(l),1)));
        A_V = -min(12*((zenithToUE - tilt_zenith)/theta_3dB)^2, SLA_V);
        A_H = -min(12*(relAzimuth/phi_3dB)^2, A_max);
        sectorGaindB(l) = G_ant_max - min(-(A_V + A_H), A_max);
    end
    %LOS state per site (SMa LOS probability, Table 7.4.2-1), shared by
    %the co-located sectors of a site
    LOSstate = false(M_sectors,1);
    invk_sum = (-log(1-r_com)*(h_com - hUT(k)) - log(1-r_res)*max(h_res - hUT(k),0))/(d_clutter*(hBS - hUT(k)));
    if r_veg > 0
        invk_sum = invk_sum - log(1-r_veg)*(h_veg - hUT(k))/(d_clutter*(hBS - hUT(k)));
    end
    for s = 1:M_sites
        d2D_site = distanceAPstoUE((s-1)*S+1); %identical for co-located sectors
        if d2D_site <= 10
            pLOS = 1;
        else
            pLOS = exp(-d2D_site*invk_sum);
        end
        LOSstate((s-1)*S+(1:S)) = rand < pLOS;
    end
    %SMa pathloss (Table 7.4.1-1) and per-link shadow fading std
    PLdB = zeros(M_sectors,1);
    sigmaSF = zeros(M_sectors,1);
    dBP = 2*pi*hBS*hUT(k)*params.fc/c0; %breakpoint distance (note 5)
    for l = 1:M_sectors
        d2D = distanceAPstoUE(l);
        d3D = distances(l,k);
        if LOSstate(l)
            if d2D < dBP
                PLdB(l) = PL1_SMa(d3D);
                sigmaSF(l) = 4;
            else
                PLdB(l) = PL1_SMa(dBP) + 40*log10(d3D/dBP);
                sigmaSF(l) = 6;
            end
        else
            PLdB(l) = 161.04 - 7.1*log10(W_street) + 7.5*log10(h_bldg) ...
                - (24.37 - 3.7*(h_bldg/hBS)^2)*log10(hBS) ...
                + (43.42 - 3.1*log10(hBS))*(log10(d3D) - 3) ...
                + 20*log10(fc_GHz) - (3.2*(log10(11.75*hUT(k)))^2 - 4.97);
            sigmaSF(l) = 8;
        end
    end
    %If this is not the first UE
    if k-1>0
        
        %Compute distances from the new prospective UE to all other UEs
        shortestDistances = zeros(k-1,1);
        
        for i = 1:k-1
            shortestDistances(i) = min(abs(UEposition - UEpositions(i) + wrapLocations));
        end
        
        %Compute conditional mean and standard deviation necessary to
        %obtain the new shadow fading realizations, when the previous
        %UEs' shadow fading realization have already been generated.
        %This computation is based on Theorem 10.2 in "Fundamentals of
        %Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn = 2.^(-shortestDistances/decorr);
        term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
        meanvalues = term1*shadowAPrealizations(1:k-1,:);
        stdvalue = sqrt(1 - term1*newcolumn);

    else %If this is the first UE

        %Add the UE and begin to store shadow fading correlation values
        meanvalues = 0;
        stdvalue = 1;
        newcolumn = [];

    end
    
    %Generate the normalized shadow fading realizations: one per site,
    %shared by the co-located sectors of that site
    shadowingSite = stdvalue*randn(1,M_sites);
    shadowing = meanvalues + repelem(shadowingSite,1,S);

    %Compute the channel gain divided by noise power: SMa pathloss,
    %sector antenna gain, per-link scaled shadow fading, and O2I loss
    gainOverNoisedB(:,k) = -PLdB + sectorGaindB + sigmaSF.*shadowing' - O2IdB(k) - noiseVariancedBm;
    %Donor-side repeater antenna: the CPE aims ONE panel at its serving
    %(donor) sector - the strongest link just computed - and every
    %sector's donor gain is that single panel's pattern response. Signals
    %from other directions are therefore suppressed by the panel pattern
    %instead of each being caught by its own best-pointing panel, so the
    %repeater does not re-amplify interference at full gain.
    if k <= K_FWA
        [~,m_donor] = max(gainOverNoisedB(:,k));
        azDonor = rad2deg(angle(APpositionsWrapped(m_donor,whichpos(m_donor)) - UEposition));
        panelAz = repOrientations(k) + (0:S-1)*360/S;
        [~,p_star] = min(abs(mod(azDonor - panelAz + 180, 360) - 180));
        boresight_p = panelAz(p_star);
        for l = 1:M_sectors
            azCPEtoBS = rad2deg(angle(APpositionsWrapped(l,whichpos(l)) - UEposition));
            zenCPEtoBS = 90 - rad2deg(atan(distanceVertical/max(distanceAPstoUE(l),1))); %pointing up at the BS
            relAzimuth = mod(azCPEtoBS - boresight_p + 180, 360) - 180;
            A_V = -min(12*((zenCPEtoBS - 90)/theta_3dB)^2, SLA_V);
            A_H = -min(12*(relAzimuth/phi_3dB)^2, A_max);
            repDonorGaindB(l,k) = G_ant_max - min(-(A_V + A_H), A_max);
        end
    end
    
    %Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:k-1,k) = newcolumn;
    shadowCorrMatrix(k,1:k-1) = newcolumn';
    shadowAPrealizations(k,:) = shadowing;
    
    %Store the UE position
    UEpositions(k) = UEposition;
    % If this is not the first CPE
    if k-1>0
        %Compute conditional mean and standard deviation necessary to
        %obtain the new shadow fading realizations, when the previous
        %UEs' shadow fading realization have already been generated.
        %This computation is based on Theorem 10.2 in "Fundamentals of
        %Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn_ue = (sigma_sf^2).*2.^(-abs(UEposition - UEpositions(1:k-1))/decorr);
        term1 = newcolumn_ue'/shadowCorrMatrix_ue(1:k-1,1:k-1);
        meanvalues_ue = term1*shadowCPErealizations_ue(1:k-1,:);
        stdvalue_ue = sqrt(sigma_sf^2 - term1*newcolumn_ue);

    else %If this is the first UE

        %Add the UE and begin to store shadow fading correlation values
        meanvalues_ue = 0;
        stdvalue_ue = sigma_sf;
        newcolumn_ue = [];

    end
    shadowing_ue = meanvalues_ue + stdvalue_ue*randn(1,K);
    %Service-side repeater antenna gain: transmitter CPE l reaches user k
    %through its best-pointing sector panel (boresight at the horizon);
    %non-CPE transmitters (rows K_FWA+1:K) stay omnidirectional
    repServiceGaindB = zeros(K,1);
    for l = 1:K_FWA
        if l ~= k
            azCPEtoUE = rad2deg(angle(UEposition - UEpositions(l)));
            d2D_ue = abs(UEposition - UEpositions(l));
            zenCPEtoUE = 90 + rad2deg(atan((hUT(l) - hUT(k))/max(d2D_ue,1)));
            bestGain = -Inf;
            for s = 1:S
                relAzimuth = mod(azCPEtoUE - (repOrientations(l) + (s-1)*360/S) + 180, 360) - 180;
                A_V = -min(12*((zenCPEtoUE - 90)/theta_3dB)^2, SLA_V);
                A_H = -min(12*(relAzimuth/phi_3dB)^2, A_max);
                bestGain = max(bestGain, G_ant_max - min(-(A_V + A_H), A_max));
            end
            repServiceGaindB(l) = bestGain;
        end
    end
    %Compute the channel gain divided by noise power; the O2I loss of an
    %indoor receiving UE also applies to its inter-UE links
    gainOverNoisedB_ue(:,k) = constantTerm - alpha*log10(distancesUEs(:,k)) + repServiceGaindB + shadowing_ue' - O2IdB(k) - noiseVariancedBm;

    % %Update shadowing correlation matrix and store realizations
    shadowCorrMatrix_ue(1:k-1,k) = newcolumn_ue;
    shadowCorrMatrix_ue(k,1:k-1) = newcolumn_ue';
    shadowCPErealizations_ue(k,:) = shadowing_ue;  

    %Inter-UE spatial correlation (i.i.d. fading): scaled identities
    for l = 1:K_FWA
        R_interue(:,:,l,k) = db2pow(gainOverNoisedB_ue(l,k))*eye(N_UE_FWA);
    end
    
    %BS-side spatial correlation (i.i.d. fading): scaled identities
    for l = 1:M_sectors
        R_gNB(:,:,l,k) = db2pow(gainOverNoisedB(l,k))*eye(N);
        if (k<=K_FWA)
            R_cpe(:,:,l,k) = eye(N_UE_FWA);
        else
            R_ue(:,:,l,k-K_FWA) = eye(N_UE_cell);
        end
    end
end



%Sector association: gainOverNoisedB includes the sector antenna pattern
%gain, so picking the strongest BS entries maps each CPE/UE to the
%sector(s) whose beam covers it
for k = 1:K_FWA
    [~, idxs] = sort(gainOverNoisedB(:,k), 'descend');
    idxs_not_chosen = idxs((Lmax+1):end);
    D_FWA(idxs_not_chosen,k) = 0;
end
for k = 1:K-K_FWA
    [~, idxs] = sort(gainOverNoisedB(:,k+K_FWA), 'descend');
    idxs_not_chosen = idxs(2:end);
    D_cell(idxs_not_chosen,k) = 0;
end
end