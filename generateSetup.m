function [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,pilotIndex,D_FWA,D_cell,APpositions,UEpositions,distances,distancesUEs,isIndoor,hUT,O2IdB,repDonorGaindB,mobState] = generateSetup(params,seed,mobState)
%Generates one realization of the suburban (SMa) simulation setup:
%sectored BS antenna gains (TR 38.901 Table 7.3-1), SMa pathloss / LOS /
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

%Temporal mobility mode: gated by params.TEMPORAL_MOBILITY (explicit,
%caller-visible flag) AND the presence of a mobState from a previous
%snapshot. When both hold, the latent large-scale variables (shadow
%fading, LOS states, O2I losses, indoor floor heights, repeater panel
%orientations) are EVOLVED consistently with each user's displacement
%instead of being redrawn, following the spatial-consistency intent of
%TR 38.901 Clause 7.6.3. Stationary users (the FWA CPEs, zero
%displacement) keep their values exactly; the caller owns the RNG stream
%across snapshots. With the flag off (the default, e.g. SimulationMain.m)
%every call behaves exactly as before this feature was added, regardless
%of what is passed as mobState.
evolveMode = isfield(params,'TEMPORAL_MOBILITY') && params.TEMPORAL_MOBILITY ...
    && (nargin > 2) && ~isempty(mobState);

%Set the seed number if it is specified other than zero
if (nargin>1)&&(seed>0)&&(~evolveMode)
    rng(seed)
end

%Communication bandwidth (Hz)
B = params.Band;

%Noise figure (in dB)
noiseFigure = params.noiseFigure;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Shadow-fading correlation distances in the horizontal plane (m),
%TR 38.901 Table 7.5-6. BS links (UMa): 37 m LOS / 50 m NLOS; inter-UE
%small-cell links (UMi): 10 m LOS / 13 m NLOS. One normalized field per
%link is drawn before the LOS state applies, so the NLOS values are used
%for both the cross-user spatial correlation and the temporal AR(1)
%evolution along the mobility trajectories
decorr = 50;    %UMa BS-side field
decorr_ue = 13; %UMi inter-UE field

%LOS-state correlation distance (m) for temporal evolution, per the
%spatial-consistency procedure of TR 38.901 Clause 7.6.3 (Table
%7.6.3.1-2, macro-scenario LOS state)
d_LOS = 50;

%Gaussian local scattering: angular standard deviation (deg) of the
%scattering cluster around the nominal azimuth, applied at the gNB
%arrays (all links) and at the NCR donor arrays (BS -> CPE links).
%Closed-form correlation (Bjornson, Hoydis, Sanguinetti, "Massive MIMO
%Networks", eq. (2.24)): R(l,m) = beta * exp(1i*pi*(l-m)*sin(phi))
%* exp(-sigma^2/2 * (pi*(l-m)*cos(phi))^2) for a half-wavelength ULA -
%a Toeplitz matrix built in O(N) and factored once by Cholesky, so the
%correlated model adds NO per-realization cost (the returned R_gNB /
%R_cpe / R_interue are CHOLESKY FACTORS, not covariances)
if isfield(params,'ASD_deg')
    ASD = params.ASD_deg*pi/180;
else
    ASD = 30*pi/180; %paper Section II value
end
locScatCol = @(Nant,phi,beta) beta*exp(1i*pi*(0:Nant-1)'*sin(phi)) ...
    .*exp(-(ASD^2/2)*(pi*(0:Nant-1)'*cos(phi)).^2);

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
%LOS probability clutter parameters (Table 7.4.2-1, SMa row): the
%suburban ZONING MIX expressed in the standard's own clutter densities -
%commercial (strip malls), residential (homes and apartment complexes),
%and natural vegetation
d_clutter = 30;
h_com = 20; h_res = 8; h_veg = 15;
r_com = 0.02; r_res = 0.18; r_veg = params.r_vegetation;
%O2I building penetration (Clause 7.4.3.1): low-loss (traditional) and
%high-loss (thermally efficient, metal-coated IRR glass) building models.
%The mix is a simulation parameter (params.high_loss_ratio); ITU-R
%M.2412 evaluations use 20% high / 80% low loss in Dense Urban eMBB and
%100% low loss in Rural eMBB
L_glass = 2 + 0.2*fc_GHz;
L_IRRglass = 23 + 0.3*fc_GHz;
L_concrete = 5 + 4*fc_GHz;
PL_tw = 5 - 10*log10(0.3*10^(-L_glass/10) + 0.7*10^(-L_concrete/10));
sigma_P = 4.4;
PL_tw_high = 5 - 10*log10(0.7*10^(-L_IRRglass/10) + 0.3*10^(-L_concrete/10));
sigma_P_high = 6.5;
if isfield(params,'high_loss_ratio')
    high_loss_ratio = params.high_loss_ratio;
else
    high_loss_ratio = 0;
end
%SMa LOS pathloss below the breakpoint (Table 7.4.1-1); PL2 above the
%breakpoint is PL1 evaluated at d_BP plus 40log10(d3D/d_BP)
PL1_SMa = @(d3D) 20*log10(40*pi*d3D*fc_GHz/3) + min(0.03*h_bldg^1.72,10)*log10(d3D) ...
                 - min(0.044*h_bldg^1.72,14.77) + 0.002*log10(h_bldg)*d3D;

%Inter-UE (CPE/NCR -> user) links: the NCR access link is modelled as a
%SMALL CELL using the TR 38.901 UMi-Street-Canyon models (Table 7.4.1-1
%pathloss, Table 7.4.2-1 LOS probability, sigma_SF 4 dB LOS / 7.82 dB
%NLOS), with the top-floor-mounted CPE in the micro-BS role (UMi assumes
%hBS = 10 m; our CPE mounts span 10.5-22.5 m and enter the formulas
%through the actual heights). This removes the last non-3GPP propagation
%model from the setup.
PL1_UMi = @(d3D) 32.4 + 21*log10(d3D) + 20*log10(fc_GHz);
PL2_UMi = @(d3D,dBP,h1,h2) 32.4 + 40*log10(d3D) + 20*log10(fc_GHz) ...
                 - 9.5*log10(dBP^2 + (h1 - h2)^2);
PLN_UMi = @(d3D,h2) 22.4 + 35.3*log10(d3D) + 21.3*log10(fc_GHz) - 0.3*(h2 - 1.5);
sigma_UMi_LOS = 4;
sigma_UMi_NLOS = 7.82;

%Prepare to save results
gainOverNoisedB = zeros(M_sectors,K);
gainOverNoisedB_ue = zeros(K,K);
R_gNB = zeros(N,N,M_sectors,K);
R_cpe = zeros(N_UE_FWA,N_UE_FWA,M_sectors,K_FWA);
%Only CPEs (indices 1:K_FWA) act as inter-UE transmitters, so the third
%dimension is K_FWA rather than K (a K x K allocation would waste GBs at
%the deployment UE counts)
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
if evolveMode
    %per-user displacement since the previous snapshot: drives the
    %temporal correlation of every latent large-scale variable
    deltaD = abs(UEpositions - mobState.UEpositions);
end
%per-site latent stores returned in mobState for the next snapshot
shadowSiteStore = zeros(K,M_sites);
LOSsiteStore = false(M_sites,K);
LOSueStore = false(K_FWA,K); %LOS state of each CPE/NCR -> user link (UMi model)
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

shadowCorrMatrix_ue = ones(K,K); %normalized field; per-link sigma applied at use
shadowCPErealizations_ue = zeros(K,K);

%Indoor state, height, and O2I penetration loss of each user. CPEs are
%outdoor-mounted at the top floor of their building. Cellular UEs are indoors with
%probability params.indoor_UT_ratio (Table 7.2-1), on a uniformly drawn
%floor of a 4-8 floor building (UMa floor model), and receive a
%UT-specific O2I loss (low/high-loss models of Clause 7.4.3.1) applied
%to all their links
isIndoor = false(K,1);
if isfield(params,'isIndoorUE') && ~isempty(params.isIndoorUE)
    %Location-determined indoor state supplied by the caller: UEs inside
    %the housing complex are indoors, UEs outside it are outdoors
    isIndoor(K_FWA+1:K) = params.isIndoorUE;
else
    %Fall back to an i.i.d. draw with the Table 7.2-5 indoor UT ratio
    isIndoor(K_FWA+1:K) = rand(K-K_FWA,1) < params.indoor_UT_ratio;
end
if evolveMode
    %indoor users stay inside the same building while walking: floor
    %height and O2I penetration loss are fixed across snapshots
    hUT = mobState.hUT;
    O2IdB = mobState.O2IdB;
else
    %Zone-dependent building heights. Apartment complexes: 2-4 floors;
    %single-family homes: 1-2 floors; strip mall: single-storey
    %commercial. Zone of each cellular UE in params.ue_zone (1 =
    %apartment, 2 = home, 3 = mall); mall-mounted CPEs flagged in
    %params.cpe_is_mall. Fallbacks keep older callers working.
    if isfield(params,'apt_floors'), apt_flr = params.apt_floors; else, apt_flr = [4 8]; end
    if isfield(params,'home_floors'), home_flr = params.home_floors; else, home_flr = [1 2]; end
    if isfield(params,'op_floors'), op_flr = params.op_floors; else, op_flr = [2 3]; end
    if isfield(params,'ue_zone') && ~isempty(params.ue_zone)
        ue_zone = params.ue_zone;
    else
        ue_zone = ones(K-K_FWA,1); %default: all apartment-zone
    end
    if isfield(params,'cpe_zone') && ~isempty(params.cpe_zone)
        cpe_zone = params.cpe_zone;
    else
        cpe_zone = ones(K_FWA,1); %default: all apartment-mounted
    end
    hUT = params.hr*ones(K,1);
    %CPE mounting, all OUTDOOR: apartment CPEs at the TOP FLOOR of their
    %2-4 floor building; single-family CPEs at the top floor of their
    %1-2 floor home (rooftop/eave mount); strip-mall CPEs on the
    %single-storey rooftop (~6 m) with open sightlines
    for k = 1:K_FWA
        switch cpe_zone(k)
            case 4 %office park: top floor of a low-rise office building
                Nfl = randi([op_flr(1) op_flr(2)]);
                hUT(k) = 3*(Nfl - 1) + 1.5;
            case 3 %strip mall
                hUT(k) = 6.0;
            case 2 %single-family home
                Nfl = randi([home_flr(1) home_flr(2)]);
                hUT(k) = 3*(Nfl - 1) + 1.5;
            otherwise %apartment
                Nfl = randi([apt_flr(1) apt_flr(2)]);
                hUT(k) = 3*(Nfl - 1) + 1.5;
        end
    end
    O2IdB = zeros(K,1);
    for k = K_FWA+1:K
        if isIndoor(k)
            d2Din = min(10*rand, 10*rand); %depth: min of two U(0,10) m
            switch ue_zone(k-K_FWA)
                case 4 %office park: low-rise, HIGH-LOSS commercial (same
                       %M.2412 commercial-glass observation as the mall)
                    Nfl = randi([op_flr(1) op_flr(2)]);
                    nfl = randi([1 Nfl]);
                    hUT(k) = 3*(nfl - 1) + 1.5;
                    O2IdB(k) = PL_tw_high + 0.5*d2Din + sigma_P_high*randn;
                case 3 %strip mall: single storey, HIGH-LOSS commercial
                       %(metal-coated IRR glass "currently appears to be
                       %more predominant in commercial buildings than in
                       %residential buildings" - ITU-R M.2412 / TR 38.901
                       %Clause 7.4.3.1)
                    hUT(k) = 1.5;
                    O2IdB(k) = PL_tw_high + 0.5*d2Din + sigma_P_high*randn;
                case 2 %single-family home: 1-2 floors, low-loss
                    Nfl = randi([home_flr(1) home_flr(2)]);
                    nfl = randi([1 Nfl]);
                    hUT(k) = 3*(nfl - 1) + 1.5;
                    O2IdB(k) = PL_tw + 0.5*d2Din + sigma_P*randn;
                otherwise %apartment complex: 2-4 floors; residential
                          %buildings are low-loss unless the caller sets
                          %a high-loss ratio
                    Nfl = randi([apt_flr(1) apt_flr(2)]);
                    nfl = randi([1 Nfl]);
                    hUT(k) = 3*(nfl - 1) + 1.5;
                    if rand < high_loss_ratio
                        O2IdB(k) = PL_tw_high + 0.5*d2Din + sigma_P_high*randn;
                    else
                        O2IdB(k) = PL_tw + 0.5*d2Din + sigma_P*randn;
                    end
            end
        end
    end
end

%Sectored antennas at the CPEs / network-controlled repeaters: every CPE
%(each one is a potential repeater) carries S sector panels with the same
%element pattern as the gNB sectors (Table 7.3-1), boresight at the
%horizon, at a random installation azimuth. A link goes through the
%best-pointing sector. Donor-side gains (CPE towards each sector BS) are
%returned in repDonorGaindB; service-side gains (CPE towards each user)
%are embedded in gainOverNoisedB_ue below.
if evolveMode
    repOrientations = mobState.repOrientations; %fixed CPE hardware installation
else
    repOrientations = 360*rand(K_FWA,1); %installation azimuth per CPE (deg)
end
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
    relAzBS = zeros(M_sectors,1); %UE bearing relative to each sector boresight (deg)
    for l = 1:M_sectors
        %azimuth of UE k relative to sector boresight, wrapped to [-180,180]
        azimuthToUE = rad2deg(angle(UEposition - APpositionsWrapped(l,whichpos(l))));
        relAzimuth = mod(azimuthToUE - boresights(l) + 180, 360) - 180;
        relAzBS(l) = relAzimuth;
        %zenith angle of the LOS direction at the BS (90 deg = horizon)
        zenithToUE = 90 + rad2deg(atan(distanceVertical/max(distanceAPstoUE(l),1)));
        A_V = -min(12*((zenithToUE - tilt_zenith)/theta_3dB)^2, SLA_V);
        A_H = -min(12*(relAzimuth/phi_3dB)^2, A_max);
        sectorGaindB(l) = G_ant_max - min(-(A_V + A_H), A_max);
    end
    %LOS state per site (SMa LOS probability, Table 7.4.2-1: exponential
    %clutter model over the commercial/residential/vegetation mix),
    %shared by the co-located sectors of a site, sticky under mobility
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
        if evolveMode && rand >= min(deltaD(k)/d_LOS,1)
            %sticky LOS state within the correlation distance d_LOS
            LOSstate((s-1)*S+(1:S)) = mobState.LOSstate(s,k);
        else
            LOSstate((s-1)*S+(1:S)) = rand < pLOS;
        end
    end
    LOSsiteStore(:,k) = LOSstate(1:S:M_sectors);
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
    if evolveMode
        %AR(1) evolution of the per-site normalized shadow fading using
        %the same spatial correlation function applied across users,
        %2^(-d/decorr): exact along each user's own trajectory (a
        %stationary user keeps its value); the cross-user correlation is
        %inherited from the initial snapshot
        a_sf = 2^(-deltaD(k)/decorr);
        shadowingSite = a_sf*mobState.shadowSite(k,:) + sqrt(1-a_sf^2)*randn(1,M_sites);
        shadowing = repelem(shadowingSite,1,S);
        newcolumn = [];
    else
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
    end
    shadowSiteStore(k,:) = shadowing(1:S:M_sectors);

    %Compute the channel gain divided by noise power: SMa pathloss,
    %sector antenna gain, per-link scaled shadow fading, and O2I loss
    gainOverNoisedB(:,k) = -PLdB + sectorGaindB + sigmaSF.*shadowing' - O2IdB(k) - noiseVariancedBm;
    %Donor-side repeater antenna: the CPE aims ONE panel at its serving
    %(donor) sector - the strongest link just computed - and every
    %sector's donor gain is that single panel's pattern response. Signals
    %from other directions are therefore suppressed by the panel pattern
    %instead of each being caught by its own best-pointing panel, so the
    %repeater does not re-amplify interference at full gain.
    relAzCPE = zeros(M_sectors,1);
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
            relAzCPE(l) = relAzimuth; %BS bearing relative to the donor panel
            A_V = -min(12*((zenCPEtoBS - 90)/theta_3dB)^2, SLA_V);
            A_H = -min(12*(relAzimuth/phi_3dB)^2, A_max);
            repDonorGaindB(l,k) = G_ant_max - min(-(A_V + A_H), A_max);
        end
    end
    
    %Update shadowing correlation matrix and store realizations
    if ~evolveMode
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
    end
    shadowAPrealizations(k,:) = shadowing;
    
    %Store the UE position
    UEpositions(k) = UEposition;
    if evolveMode
        %per-link AR(1) evolution of the inter-UE shadow fading: both
        %endpoints' displacements decorrelate a link (the transmitting
        %CPEs are stationary, so only the receiver's motion matters)
        a_ue = 2.^(-(deltaD(k) + deltaD)'/decorr_ue); %1 x K
        shadowing_ue = a_ue.*mobState.shadow_ue(k,:) + sqrt(1-a_ue.^2).*randn(1,K);
        newcolumn_ue = [];
    else
    % If this is not the first CPE
    if k-1>0
        %Compute conditional mean and standard deviation necessary to
        %obtain the new shadow fading realizations, when the previous
        %UEs' shadow fading realization have already been generated.
        %This computation is based on Theorem 10.2 in "Fundamentals of
        %Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn_ue = 2.^(-abs(UEposition - UEpositions(1:k-1))/decorr_ue);
        term1 = newcolumn_ue'/shadowCorrMatrix_ue(1:k-1,1:k-1);
        meanvalues_ue = term1*shadowCPErealizations_ue(1:k-1,:);
        stdvalue_ue = sqrt(1 - term1*newcolumn_ue);

    else %If this is the first UE

        %Add the UE and begin to store shadow fading correlation values
        meanvalues_ue = 0;
        stdvalue_ue = 1;
        newcolumn_ue = [];

    end
    shadowing_ue = meanvalues_ue + stdvalue_ue*randn(1,K);
    end
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
    %UMi-Street-Canyon pathloss and LOS state per CPE/NCR -> user link
    %(the NCR access link as a small cell). Rows K_FWA+1:K (UE -> UE) are
    %never consumed downstream and use the NLOS expression directly
    PL_ue = zeros(K,1);
    sigmaSF_ue = sigma_UMi_NLOS*ones(K,1);
    for l = 1:K
        d2D_ue = max(distancesUEs(l,k), 1);
        d3D_ue = sqrt(d2D_ue^2 + (hUT(l) - hUT(k))^2);
        dBP_ue = 4*(hUT(l) - 1)*(hUT(k) - 1)*params.fc/c0; %d'BP, hE = 1 (UMi)
        if d2D_ue < dBP_ue
            PL_LOS_ue = PL1_UMi(d3D_ue);
        else
            PL_LOS_ue = PL2_UMi(d3D_ue,dBP_ue,hUT(l),hUT(k));
        end
        if l <= K_FWA && l ~= k
            %LOS probability (Table 7.4.2-1, UMi row) with a sticky LOS
            %state along the receiver's trajectory (the CPE is static)
            if d2D_ue <= 18
                pLOS_ue = 1;
            else
                pLOS_ue = 18/d2D_ue + exp(-d2D_ue/36)*(1 - 18/d2D_ue);
            end
            if evolveMode && rand >= min(deltaD(k)/d_LOS,1)
                LOS_l = mobState.LOS_ue(l,k);
            else
                LOS_l = rand < pLOS_ue;
            end
            LOSueStore(l,k) = LOS_l;
            if LOS_l
                PL_ue(l) = PL_LOS_ue;
                sigmaSF_ue(l) = sigma_UMi_LOS;
            else
                PL_ue(l) = max(PL_LOS_ue, PLN_UMi(d3D_ue,hUT(k)));
            end
        else
            PL_ue(l) = max(PL_LOS_ue, PLN_UMi(d3D_ue,hUT(k)));
        end
    end
    %Compute the channel gain divided by noise power; the O2I loss of an
    %indoor receiving UE also applies to its inter-UE links
    gainOverNoisedB_ue(:,k) = -PL_ue + repServiceGaindB + sigmaSF_ue.*shadowing_ue' - O2IdB(k) - noiseVariancedBm;

    % %Update shadowing correlation matrix and store realizations
    if ~evolveMode
        shadowCorrMatrix_ue(1:k-1,k) = newcolumn_ue;
        shadowCorrMatrix_ue(k,1:k-1) = newcolumn_ue';
    end
    shadowCPErealizations_ue(k,:) = shadowing_ue;  

    %Inter-UE links keep i.i.d. fading (handset/CPE service arrays in
    %rich local clutter): the factor of beta*I is sqrt(beta)*I
    for l = 1:K_FWA
        R_interue(:,:,l,k) = sqrt(db2pow(gainOverNoisedB_ue(l,k)))*eye(N_UE_FWA);
    end
    
    %BS-side spatial correlation: Gaussian local scattering around the
    %link's nominal azimuth at the gNB ULA (all links) and at the NCR
    %donor array (BS -> CPE links); stored as CHOLESKY FACTORS. Cellular
    %handsets (R_ue) stay uncorrelated (rich local scattering).
    for l = 1:M_sectors
        beta_lk = db2pow(gainOverNoisedB(l,k));
        Rcol = locScatCol(N, deg2rad(relAzBS(l)), beta_lk);
        Rfull = toeplitz(Rcol, Rcol');
        R_gNB(:,:,l,k) = chol(Rfull + 1e-9*beta_lk*eye(N), 'lower');
        if (k<=K_FWA)
            Rcol = locScatCol(N_UE_FWA, deg2rad(relAzCPE(l)), 1);
            Rfull = toeplitz(Rcol, Rcol');
            R_cpe(:,:,l,k) = chol(Rfull + 1e-9*eye(N_UE_FWA), 'lower');
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

%Latent-variable state handed to the next temporal-mobility snapshot
mobState = struct('UEpositions',UEpositions,'shadowSite',shadowSiteStore, ...
    'LOSstate',LOSsiteStore,'LOS_ue',LOSueStore,'O2IdB',O2IdB,'hUT',hUT, ...
    'repOrientations',repOrientations,'shadow_ue',shadowCPErealizations_ue);
end