close all;
clear;
tStart = tic;
aID = getenv('SLURM_ARRAY_TASK_ID');

% This is for running on a cluster in parallel
% the bash script should give the aID as input
if (isempty(aID))
    warning('aID is empty. Trying SLURM ID.')
    aID = getenv('SLURM_ARRAY_TASK_ID');
end
if(isempty(aID))
    warning('aID is empty. Replacing it with 0010.')
    aID = '0084';
end
%RNG seed.
rng(str2double(aID),'twister');
%% Mobility params (TR 38.901 Table 7.2-5)
params.MOBILE = 1;
params.ue_velocity_outdoor = 40/3.6; %outdoor (in-car) UTs: 40 km/h
params.ue_velocity_indoor = 3/3.6;   %indoor UTs: 3 km/h
params.Ts = 1.67e-5; %channel sample time (s)
params.fc = 3.7e9; %carrier frequency (Hz)
params.c = 3e8; %speed of light (m/s)
%params.mob_rho (one channel-aging factor per cellular UE) is set after
%generateSetup, once the indoor/outdoor state of each UE is known
%% Rate requirements
%FCC 2024 Section 706 Report (FCC 24-27, https://docs.fcc.gov/public/attachments/FCC-24-27A1.pdf):
% - mobile 5G-NR coverage benchmark 35/3 Mbps DL/UL
% - fixed broadband benchmark 100/20 Mbps
% - long-term fixed broadband goal 1000/500 Mbps
%Ericsson Mobility Report (FWA): busy-hour sustained per active user ~5 Mbps
% (https://www.ericsson.com/en/reports-and-papers/mobility-report/articles/fixed-wireless-access)
params.r_min_cell = 10e6; %busy-hour sustained cell-UE floor (~5-10 Mbps; the FCC 35/3
                          %figure is a coverage benchmark, not a per-user busy-hour target)
params.r_max_FWA = 1e9;   %FWA rate cap = FCC long-term fixed broadband goal (1 Gbps),
                          %also the Verizon 5G Home "Ultimate" gigabit tier
%%
params.BEAM = 0;
%%


numCPE_all = 50; %CPEs per site
Band = 100e6;
params.noiseFigure = 9; %receiver noise figure (dB)

%% Suburban Macro (SMa) deployment, 3GPP TR 38.901 Table 7.2-5
params.ISD = 1299; %inter-site distance (m); 1732 m also listed, ISDs of 1200-1800 m valid (note 2)
params.deployRange = params.ISD/2; %cell radius (half ISD); sets the wrap-around span
params.complexRadius = 200; %housing complex radius around each site (m): CPEs and the indoor
                            %UEs live here; outdoor (in-car) UEs roam the rest of the cell
params.coverageRange_sub6 = 2*params.deployRange + 30; %wrap-around square length
params.min_dist_2D = 35; %minimum BS-UT 2D distance (m)
params.num_antennas_per_gNB = 64;
params.rho_tot = 10^(0.1*49); %BS Tx power 49 dBm per sector, per the SMa calibration
                              %assumptions of TR 38.901 clause 7.8 (ETSI TR 138 901 V19.2.0)
params.CELL_REPEAT = 1;
%Number of antennas per UE
params.N_UE_FWA = 8;
params.N_UE_cell = 2; %4;
params.hr = 1.5; %outdoor UT height (m)
params.hr_cpe = 4.5; %CPE at top floor of 2-storey residential building (m)
params.ht_bs = 35; %BS antenna height (m)
%% Sectored BS antennas: 3 sectors per site, element pattern per TR 38.901 Table 7.3-1
params.sectors_per_site = 3;
sector_offsets = [30; 150; 270]; %sector boresight azimuths within a site (deg)
params.theta_3dB = 65; %vertical 3 dB beamwidth (deg)
params.phi_3dB = 65; %horizontal 3 dB beamwidth (deg)
params.SLA_V = 30; %side-lobe attenuation, vertical (dB)
params.A_max = 30; %front-back attenuation (dB)
params.G_ant_max = 8; %max element gain (dBi)
%Mechanical downtilt in zenith-angle convention (90 deg = horizon), per
%TR 38.901 SMa calibration: 92 deg for ISD = 1732 m, 95 deg for ISD = 1299 m
if params.ISD <= 1500
    params.tilt_zenith = 95;
else
    params.tilt_zenith = 92;
end
%% SMa propagation, TR 38.901 Clause 7.4 (pathloss, LOS probability, O2I)
params.h_bldg = 10; %average building height h (m), SMa default
params.W_street = 10; %average street width W (m), SMa default
params.r_vegetation = 0.10; %vegetation clutter density: 0 none, 0.10 sparse, 0.20 dense
params.indoor_UT_ratio = 0.8; %fraction of cellular UEs indoors (Table 7.2-5)
lambda_BS = 5;
lambda_SC = 0; %no small cells (kept for the result-file format)
%Active-user densities (busy-hour CONCURRENTLY ACTIVE UEs, not
%subscribers). Suburban residential density averages ~1050 people/km^2
%(Demographia international suburban densities, http://demographia.com/db-intlsub.htm;
%FiveThirtyEight, https://fivethirtyeight.com/features/how-suburban-are-big-american-cities/),
%rising to ~2000-4000/km^2 in a housing complex; with a ~3-5% busy-hour
%active share this gives ~80-200 active UEs/km^2 indoors. Suburban roads
%carry ~4-15 active in-car users/km^2. The resulting ~9 active UEs per
%sector matches the full-buffer evaluation convention of 10 UEs per TRxP
%(ITU-R M.2412-0, https://www.itu.int/dms_pub/itu-r/opb/rep/R-REP-M.2412-2017-PDF-E.pdf;
%3GPP TR 38.913 / ETSI TR 138 913), and the 150/4 defaults reproduce the
%80:20 indoor:outdoor UT ratio of TR 38.901 Table 7.2-5.
lambda_UE_res_arr = 150; %50:50:400; %active UEs per km^2, residential block
lambda_UE_road_arr = 4;  %2:2:20;    %active in-car UEs per km^2, roads outside (paired)
params.Lmax = 1;
params.preLogFactor = 1;
params.loss_pc_cell = 5/100; %enforce r_min_cell at the 5th-percentile cell UE
SI_cancel_arr = 0; %-20:5:0; %SI cancellation factor sweep (dB) for the FWA phase;
                   %set to -20:5:0 to reproduce the retired SE_comparison.m experiment
params.HW_IMPAIRMENTS = 1;  % 1 = hardware impairments on, 0 = ideal hardware
%EVM-based impairment factors: an error vector magnitude of e turns a
%fraction e^2 of the signal power into distortion, so K = 1 - e^2.
%3GPP minimum signal quality: BS transmit EVM 3.5% for 256QAM / 8% for
%64QAM, TS 38.104 Table 6.5.2-1 (https://www.etsi.org/deliver/etsi_ts/138100_138199/138104/);
%UE EVM uses the same classes (TS 38.101-1); NR repeater conformance
%reuses them too (TS 38.115-2, https://www.etsi.org/deliver/etsi_ts/138100_138199/13811502/17.03.00_60/ts_13811502v170300p.pdf;
%RF requirements in TS 38.106). 256QAM-grade (3.5%) is assumed for the
%BS transmit chain, 64QAM-grade (8%) for the UE receive chain and for
%each half of the repeater amplify-and-forward chain.
params.Kt = 1 - 0.035^2;    % BS transmit chain, 256QAM EVM 3.5% (TS 38.104)
params.Kr = 1 - 0.08^2;     % UE receive chain, 64QAM-grade EVM 8%
params.Kt_rep = 1 - 0.08^2; % repeater transmit chain, EVM 8% (TS 38.115-2)
params.Kr_rep = 1 - 0.08^2; % repeater receive chain, EVM 8% (TS 38.115-2)
params.IMPERFECT_CSI = 1;   % 1 = beamforming uses MMSE channel estimates from noisy
                            %     CSI-RS observations, 0 = perfect CSI
params.csi_rs_offset_dB = 0;% CSI-RS EPRE relative to PDSCH EPRE (powerControlOffset,
                            % TS 38.214; standardised range -8..15 dB, default 0)
%Network-controlled repeater (NCR) amplification, per 3GPP: fixed
%amplification gain of 90 dB subject to the maximum output power
%constraint of 33 dBm (evaluation assumptions of TR 38.867, "Study on NR
%network-controlled repeaters"; RF repeater requirements in TS 38.106).
%The cap is applied per repeater at its actual operating point below, so
%the effective gain is min(G_max, P_max/P_in)
rep_gain_arr = 90; %NCR amplification gain G_max (dB), TR 38.867
params.rep_max_pow_dBm = 33; %NCR maximum output power (dBm), TR 38.867 / TS 38.106
num_rep_arr = 1:1:6; %total enabled repeaters; attachment is per donor sector,
                     %so coverage saturates at one repeater per sector (6)
params.num_repeater_per_cpe = 1; %repeaters associated per served user (fixed, not swept)
params.rep_assist_frac = 0.2;    %fraction of each sector's WEAKEST cellular UEs (cell edge,
                                 %by serving-link large-scale gain) eligible for repeater
                                 %assistance; the rest keep the pure direct path
%Number of channel realizations per setup
params.nbrOfRealizations = 10;

%Minimum guaranteed FWA rate per subscriber-plan tier (swept). Values
%span the US 5G-FWA plan landscape (Q4 2024): 50 Mbps entry floor,
%100 Mbps = FCC fixed-broadband benchmark and common base plan,
%150-200 Mbps = current median delivered speeds (Verizon/AT&T ~150,
%T-Mobile ~205), 300 Mbps = Verizon 5G Home standard-plan soft cap;
%the gigabit "Ultimate" tier is captured by params.r_max_FWA = 1 Gbps.
%Refs: Light Reading (https://www.lightreading.com/fixed-wireless-access/verizon-and-t-mobile-fwa-speeds-are-actually-increasing),
%tecknexus (https://tecknexus.com/t-mobile-verizon-and-att-boost-5g-fwa-speeds-and-subscribers/),
%Verizon FWA guide (https://www.verizon.com/home/internet/guides/what-is-fixed-wireless-access-fwa-technology/).
r_min_arr = 1e6*(50:50:300);
%% Simulation FR1 setup
for idxBSDensity = 1:length(lambda_BS)
    %% gNB locations: macro sites, each split into 3 co-located sector BSs
    params.M_sites = 2; %number of macro sites
    M_sites = params.M_sites;
    S = params.sectors_per_site;
    params.numGNB = M_sites*S; %each sector is an independent BS entry
    M_sectors = params.numGNB;
    Rsite = [0; params.ISD]; %site distances from origin (adjacent sites at ISD)
    angleSite = 2*pi*rand(M_sites,1);
    siteLocations = [Rsite.*cos(angleSite), Rsite.*sin(angleSite)];
    params.siteLocations = siteLocations;
    params.locationsBS = kron(siteLocations, ones(S,1)); %sectors co-located at their site
    params.sector_boresights = repmat(sector_offsets, M_sites, 1); %boresight azimuth per BS entry (deg)
    %% CPE locations: numCPE_all per site, uniform in the housing complex
    %annulus [min_dist_2D, complexRadius] around each site
    numCPE_tot = M_sites*numCPE_all;
    RCPE = sqrt(params.min_dist_2D^2 + (params.complexRadius^2 - params.min_dist_2D^2)*rand(numCPE_tot,1));
    angleCPE = 2*pi*rand(numCPE_tot,1);
    CPE_locations = kron(siteLocations, ones(numCPE_all,1)) + [RCPE.*cos(angleCPE), RCPE.*sin(angleCPE)];
    for idxUEDensity = 1:length(lambda_UE_res_arr)
        %% UE locations: two-density model. Active residential UEs lie in
        %the housing complex annulus [min_dist_2D, complexRadius] and are
        %indoors; active in-car UEs lie on the roads in the outer ring
        %[complexRadius, deployRange] and are outdoors. Counts follow the
        %busy-hour active-user densities applied to the actual areas
        lambda_road = lambda_UE_road_arr(min(idxUEDensity,end));
        A_complex = pi*((params.complexRadius/1000)^2 - (params.min_dist_2D/1000)^2); %residential annulus (km^2)
        A_road = pi*((params.deployRange/1000)^2 - (params.complexRadius/1000)^2); %outer ring (km^2)
        numUE_in = ceil(lambda_UE_res_arr(idxUEDensity)*A_complex/S); %per sector, residential block (indoor)
        numUE_out = ceil(lambda_road*A_road/S); %per sector, roads (in-car)
        params.numUE = numUE_in + numUE_out;
        numUE = params.numUE;
        UE_locations = zeros(M_sectors*numUE,2);
        isIndoorUE = false(M_sectors*numUE,1);
        for idxSite = 1:M_sites
            Nin = S*numUE_in;
            Nout = S*numUE_out;
            Rin = sqrt(params.min_dist_2D^2 + (params.complexRadius^2 - params.min_dist_2D^2)*rand(Nin,1));
            Rout = sqrt(params.complexRadius^2 + (params.deployRange^2 - params.complexRadius^2)*rand(Nout,1));
            angleUE = 2*pi*rand(Nin+Nout,1);
            RUE = [Rin; Rout];
            base = (idxSite-1)*(Nin+Nout);
            UE_locations(base+(1:Nin+Nout),:) = siteLocations(idxSite,:) + [RUE.*cos(angleUE), RUE.*sin(angleUE)];
            isIndoorUE(base+(1:Nin)) = true;
        end
        params.UE_locations = UE_locations;
        params.isIndoorUE = isIndoorUE;
        params.numCPE = numCPE_tot;
        params.CPE_locations = CPE_locations;
        params.Band = Band; %Communication bandwidth
        [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,pilotIndex,D_FWA,D_cell,APpositions,UEpositions,distances,distancesCPEs,isIndoor,hUT,O2IdB,repDonorGaindB] = generateSetup(params,str2double(aID));
        params.rep_donor_gain = db2pow(repDonorGaindB); %donor-side repeater sector antenna gain (linear power, M_sectors x K_FWA)
        params.gainOverNoise_lin = db2pow(gainOverNoisedB); %BS-to-user large-scale gain (linear), single source of truth
        %Total downlink power arriving at each repeater's donor panel (mW),
        %summed over all sector BSs: sets the operating point of the NCR
        %output-power cap. gainOverNoise_lin is normalized by the full-band
        %noise power, which is re-applied here
        noise_mW = 10^(0.1*(-174 + 10*log10(Band) + params.noiseFigure));
        P_in_rep_mW = params.rho_tot*noise_mW*sum(params.gainOverNoise_lin(:,1:numCPE_tot).*params.rep_donor_gain,1)';
        %Per-UE channel aging factor: indoor UEs move at 3 km/h, outdoor
        %UEs at 40 km/h (TR 38.901 Table 7.2-5)
        v_ue = params.ue_velocity_outdoor*ones(M_sectors*numUE,1);
        v_ue(isIndoor(numCPE_tot+1:end)) = params.ue_velocity_indoor;
        params.mob_rho = besselj(0,2*pi*params.Ts*v_ue*params.fc/params.c);
        params.BETA_interUE = db2pow(gainOverNoisedB_ue);
        for idxnumrep = 1:length(num_rep_arr)  
            for idxrepgain = 1:length(rep_gain_arr)
                params.repeat_gain = rep_gain_arr(idxrepgain); %G_max (dB), recorded in the CSV
                %Effective per-repeater AMPLITUDE gain: fixed gain capped by
                %the NCR maximum output power at each repeater's operating
                %point (TR 38.867: 90 dB / 33 dBm)
                params.repeat_gain_amp = sqrt(min(10^(0.1*rep_gain_arr(idxrepgain)), ...
                    10^(0.1*params.rep_max_pow_dBm)./P_in_rep_mW));
                params.num_repeater_tot = num_rep_arr(idxnumrep); %total repeaters enabled
                %Enabled repeater pool for the cellular-UE phase, chosen
                %to maximise coverage: every UE votes for the repeater it
                %would associate with (same donor x service metric as the
                %rate computation, through its serving sector), and the
                %num_repeater_tot most-voted CPEs are enabled
                BETA_FWA_assoc = db2pow(gainOverNoisedB(:,1:numCPE_tot)).*D_FWA;
                votes = zeros(numCPE_tot,1);
                for kUE = 1:M_sectors*numUE
                    m_serv = find(D_cell(:,kUE)==1,1);
                    v_rep = BETA_FWA_assoc(m_serv,:)'.*params.rep_donor_gain(m_serv,:)'.*params.BETA_interUE(1:numCPE_tot,numCPE_tot+kUE);
                    [v_best,best_rep] = max(v_rep);
                    if isfinite(v_best) && v_best > 0
                        votes(best_rep) = votes(best_rep) + 1;
                    end
                end
                [~,vote_order] = sort(votes,'descend');
                params.set_repeat = vote_order(1:min(params.num_repeater_tot,numCPE_tot));
                Band = params.Band; %Communication bandwidth
                K_FWA = params.numCPE;
                params.numUE = numUE;
                K = params.numCPE + M_sectors*params.numUE; 
                if params.CELL_REPEAT
                    params.D_FWA = D_FWA;
                    params.D_cell = D_cell;
                    params.R_gNB = R_gNB;
                    params.R_cpe = R_cpe;
                    params.R_interue = R_interue;
                    params.R_ue = R_ue; 
                    nbrOfRealizations = params.nbrOfRealizations;
                    rate_dl = zeros(K-K_FWA,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est);         
                    end
                else
                    params.numCPE = 0;
                    params.CPE_locations = [];
                    params.D = D_cell;
                    params.R_gNB = R_gNB(:,:,:,1+numCPE_tot:end);
                    params.R_cpe = [];
                    params.R_interue = [];
                    params.R_ue = R_ue; 
                    nbrOfRealizations = params.nbrOfRealizations;
                    rate_dl = zeros(K-K_FWA,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_OFDM(params, channel_dl, channel_est_dl);         
                    end
                end
                mean_rate_dl_cell = mean(rate_dl,2);
                if (quantile(mean_rate_dl_cell,params.loss_pc_cell)>=params.r_min_cell)
                    params.Band = max(params.Band*(1-params.r_min_cell/quantile(mean_rate_dl_cell,params.loss_pc_cell)),0);
                    if (params.Band > 0)
                        mean_rate_dl_cell = mean_rate_dl_cell*params.r_min_cell/quantile(mean_rate_dl_cell,params.loss_pc_cell); 
                    end
                else
                    params.Band = 0;
                end
                Band_after_cell = params.Band; %band left after serving the cellular UEs
                %SI-cancellation sweep for the FWA phase (single 0 dB value
                %by default; the retired SE_comparison.m experiment swept it)
                for idxSI = 1:length(SI_cancel_arr)
                    params.SI_cancel_factor = 10^(0.1*SI_cancel_arr(idxSI));
                    params.Band = Band_after_cell;
                    params.numUE = 0;
                    params.numCPE = numCPE_tot;
                    K_FWA = params.numCPE;
                    K = params.numCPE + M_sectors*params.numUE;
                    params.CPE_locations = CPE_locations;
                    params.D_cell = [];
                    params.D_FWA = D_FWA;
                    params.D = D_FWA;
                    CPE_idxs = 1:numCPE_tot;
                    params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
                    params.R_cpe = R_cpe;
                    params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
                    params.R_ue = []; 
                    params.set_repeat = [];
                    rate_dl = zeros(K,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    save_old_rate = rate_dl;
                    save_old_mean_FWA = mean_rate_dl_FWA;
                    Band_FWA = params.Band;
                    for idxrmin = 1:length(r_min_arr)
                        params.r_min_FWA = r_min_arr(idxrmin);
                        K_FWA_max = 0;
                        params.Band = Band_FWA;
                        [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                        K_FWA_init = sum(FWA_util>0);
                        if (params.Band > 0)
                            while any(mean_rate_dl_FWA > params.r_max_FWA)
                                CPE_idxs = find(mean_rate_dl_FWA > params.r_max_FWA);
                                params.Band = params.Band*(1-(params.r_min_FWA/min(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA))));
                                rate_dl(CPE_idxs,:) = rate_dl(CPE_idxs,:)*params.r_min_FWA/min(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA));
                                mean_rate_dl_FWA(CPE_idxs) = mean(rate_dl(CPE_idxs,:),2);
                                params.set_repeat = [params.set_repeat; CPE_idxs];
                                not_set_repeat = setdiff(1:numCPE_tot,params.set_repeat);
                                for n = 1:nbrOfRealizations
                                    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, ~, ~] = computePhysicalChannels_sub6_MIMO(params);
                                    rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                                end
                                mean_rate_dl_FWA(not_set_repeat) = mean(rate_dl(not_set_repeat,:),2);
                            end
                            [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                            K_FWA_max = sum(FWA_util>0);
                        end
                        sum_FWA_rate = sum(mean_rate_dl_FWA);
                        params.CPE_locations = CPE_locations;
                        params.numCPE = numCPE_tot;
                        params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
                        params.R_cpe = R_cpe;
                        params.R_interue = R_interue;
                        params.R_ue = []; 
                        params.D_FWA = D_FWA;
                        params.set_repeat = [];
                        params.Band = Band;
                        mean_rate_dl_FWA = save_old_mean_FWA;
                        rate_dl = save_old_rate;
                        if (Band_FWA > 0) && (K_FWA_max > 0)
                            FWA_se_out = sum_FWA_rate/(Band_FWA*K_FWA_max);
                        else
                            FWA_se_out = 0; %no band or no served CPEs: report 0 instead of 0/0
                        end

                    %% Recording the Results
                    
                        %Taking care of folder directory creation etc
                        dataFolder = 'resultData';
                        rateFolder = strcat(dataFolder,'/FWA_multi_cell_repeater_fix_comp_alloc');
                        if not(isfolder(dataFolder))
                            mkdir(dataFolder)
                        end
                        if not(isfolder(rateFolder))
                            mkdir(rateFolder)
                        end
                
                        deployRange = params.deployRange;
                        numBS = size(params.locationsBS,1);
                        result_string = strcat('/results_numFWA_',num2str(params.Lmax), 'Lmax_', num2str(numCPE_all),...
                            'CPE_',num2str(lambda_BS(idxBSDensity)),...
                            'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
                            'lambdaSC_',num2str(lambda_UE_res_arr(idxUEDensity)),...
                            'lambdaUEres_',num2str(lambda_road),'lambdaUEroad_', ...
                            num2str(deployRange),'deployRange_', ...
                            num2str(params.r_min_FWA/10^6),'min_FWA_rate', ...
                            num2str(params.repeat_gain),'repeater_gain', ...
                            num2str(params.num_repeater_tot),'num_repeater_', ...
                            num2str(SI_cancel_arr(idxSI)),'SI_cancel_', ...
                            aID);
                        recording_text_file_string = strcat(rateFolder,result_string,'.csv');
                        fileID = fopen(recording_text_file_string,'w');
                        output_categories = ['lambdaBS,','lambdaSC,','numCPE,','lambdaUE,',...
                        'lambdaUE_road,','deployRange,','r_min_cell,','r_min_FWA,','num_rep,','rep_gain,','SI_cancel,','init_FWA,','max_FWA,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
                        fprintf(fileID,output_categories);
                        formatSpec = '%d,%d,%d,%d,%d,%d,%f,%f,%d,%d,%d,%d,%d,%f,%f,%f,%f\n';
                        fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_SC(idxBSDensity),numCPE_all, ...
                        lambda_UE_res_arr(idxUEDensity),lambda_road,deployRange,params.r_min_cell,params.r_min_FWA,params.num_repeater_tot,params.repeat_gain,SI_cancel_arr(idxSI),K_FWA_init,K_FWA_max,Band,Band_FWA, sum(mean_rate_dl_cell)/(Band - Band_FWA), FWA_se_out);
                        fclose(fileID);
                        % 'lambdaUE_road,','deployRange,','r_min_cell,','r_min_FWA,','num_rep,','rep_gain,','init_FWA,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
                        % lambda_UE_res_arr(idxUEDensity),lambda_road,deployRange,params.r_min_cell,params.r_min_FWA,params.num_repeater_tot,params.repeat_gain,K_FWA_init,Band,Band_FWA, sum(mean_rate_dl_cell)/(Band - Band_FWA), sum_FWA_rate/(Band_FWA*K_FWA_init));
                    end
                end %idxSI (SI-cancellation sweep)
            end
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd) 