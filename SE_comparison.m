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
% for aID = 0:99
%RNG seed.
rng(str2double(aID),'twister');
% rng(aID,'twister');
%% GUE channel parameters
params.K_Factor = 9;         %dB -- %rician factor Ground UE  % if beta_gains=1
params.RAYLEIGH=0;   %1= rayleigh, % 0=rician
params.Perf_CSI =1;
params.cov_area = 1; %0.25; % 4; %km
%% Mobility params (TR 38.901 Table 7.2-5)
params.MOBILE = 1;
params.ue_velocity_outdoor = 40/3.6; %outdoor (in-car) UTs: 40 km/h
params.ue_velocity_indoor = 3/3.6;   %indoor UTs: 3 km/h
params.Ts = 1.67e-5; %10us
params.fc = 3.7e9; %7.9 GHz FR3
params.c = 3e8; %speed of light
params.coherence_time = 1e-3;
params.n_samples = floor(params.coherence_time/params.Ts);
%params.mob_rho (one channel-aging factor per cellular UE) is set after
%generateSetup, once the indoor/outdoor state of each UE is known
%% Economics params
params.price_cell_per_bit = 4;
params.price_FWA_per_bit = 0.14;
params.r_min_cell = 35e6;
params.r_max_FWA = 1e9;
%%
params.CH_estimation = 0;  % 1= have channel estimatio
%%
params.ASD_VALUE = 0;%[0,0.25,0.5,0.75,1];  % [0,30,10]; %
params.ASD_CORR = 0;
params.Kt_Kr_vsUE  = 0; %0.175^2; %0.175^2; %[1,2,3,4];  %to save 1=AP 0.1,UE=0.1;  2=AP 0.1,UE=0.3;  3=AP 0.3,UE=0.1

params.pilot_pow = 100;  % 0.1W   % UL pilot. power (W)
params.noiseFigure = 9; % gue
params.sigma_sf =4;    
params.SC = 1; %0 - CF, 1 - SC
%% Define simulation setup

%Angular standard deviation in the local scattering model (in radians)
params.ASD_varphi = deg2rad(30); %azimuth angle
params.ASD_theta = 0; %deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
params.p = 100;
params.UE_split = 0; %fraction in I_band
params.BEAM = 0;

numCPE_all = 50; %5:5:20;
Band = 100e6;

%Prepare to save simulation results
%% Suburban Macro (SMa) deployment, 3GPP TR 38.901 Table 7.2-5
params.ISD = 1299; %inter-site distance (m); 1732 m also listed, ISDs of 1200-1800 m valid (note 2)
params.deployRange = params.ISD/2; %cell radius (half ISD); sets the wrap-around span
params.complexRadius = 200; %housing complex radius around each site (m): CPEs and the indoor
                            %UEs live here; outdoor (in-car) UEs roam the rest of the cell
params.coverageRange_sub6 = 2*params.deployRange + 30; %wrap-around square length
params.min_dist_2D = 35; %minimum BS-UT 2D distance (m)
params.num_antennas_per_gNB = 64;
params.num_antennas_per_sc = 16;
params.rho_tot = 10^(0.1*75)*(Band/1e8);
params.rho_tot_sc = 10^(0.1*55);
%Number of antennas per UE
params.N_UE_FWA = 8;
params.N_UE_cell = 2; %4;
params.hr = 1.5; %outdoor UT height (m)
params.hr_cpe = 4.5; %CPE at top floor of 2-storey residential building (m)
params.ht_bs = 35; %BS antenna height (m)
params.ht_sc = 5;
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
lambda_BS = 5; %([5 6 7 8 9 10]).^2;
lambda_SC = 0; %([5 6 7 8 9 10]).^2;
lambda_UE = 200:200:1000; %combined UE density (per km^2) over the entire cell area
params.Lmax = 1;
params.preLogFactor = 1;
params.loss_pc_cell = 5/100;
params.loss_pc_FWA = 5/100;
%Number of channel realizations per setup
params.nbrOfRealizations = 10;
%repeater config
params.repeat_gain = 10^(0.1*(6.5+20*log10(params.fc/1e6)));
params.set_repeat = [];
params.num_repeater_per_cpe = 2;
params.CELL_REPEAT = 0;
params.FWA_REPEAT = 0;
params.HW_IMPAIRMENTS = 1;  % 1 = hardware impairments on, 0 = ideal hardware
params.Kt = 0.9;            % transmitter impairment factor (1 = ideal)
params.Kr = 0.9;            % receiver impairment factor (1 = ideal)
SI_cancel_arr = -20:5:0; %SI cancel factor in dB
%% UE angular coverage range (full 360 coverage for now)
lookAngleCell{1} = [0,360];
r_min_arr = 1e6*(25:25:300);
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
    for idxUEDensity = 1:length(lambda_UE)
        %% UE locations: two-density model over the entire cell area.
        %A fraction indoor_UT_ratio (80%, Table 7.2-5) of the UEs lies in
        %the housing complex annulus [min_dist_2D, complexRadius] and is
        %indoors; the rest lies in the outer ring [complexRadius,
        %deployRange] and is outdoors (in-car). The combined density over
        %the cell is lambda_UE per km^2
        A_cell = pi*(params.deployRange/1000)^2; %cell area (km^2)
        numUE_in = ceil(params.indoor_UT_ratio*lambda_UE(idxUEDensity)*A_cell/3); %per sector, in complex
        numUE_out = ceil((1-params.indoor_UT_ratio)*lambda_UE(idxUEDensity)*A_cell/3); %per sector, outside
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
        [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,pilotIndex,D_FWA,D_cell,APpositions,UEpositions,distances,distancesCPEs,isIndoor,hUT,O2IdB] = generateSetup(params,str2double(aID));
        %Per-UE channel aging factor: indoor UEs move at 3 km/h, outdoor
        %UEs at 40 km/h (TR 38.901 Table 7.2-5)
        v_ue = params.ue_velocity_outdoor*ones(M_sectors*numUE,1);
        v_ue(isIndoor(numCPE_tot+1:end)) = params.ue_velocity_indoor;
        params.mob_rho = besselj(0,2*pi*params.Ts*v_ue*params.fc/params.c);
        ASD_VALUE = params.ASD_VALUE;
        ASD_CORR = params.ASD_CORR;
        Kt_Kr_vsUE = params.Kt_Kr_vsUE;
        K_Factor = params.K_Factor;
        RAYLEIGH=params.RAYLEIGH;   %1= rayleigh, % 0=rician
        Perf_CSI = params.Perf_CSI;
        cov_area = params.cov_area;
        pilot_pow = params.pilot_pow; 
        noiseFigure = params.noiseFigure;
        sigma_sf = params.sigma_sf;
        Band = params.Band; %Communication bandwidth
        params.numCPE = 0;
        params.CPE_locations = [];
        K_FWA = params.numCPE;
        K = params.numCPE + M_sectors*params.numUE; 
        params.BETA = db2pow(gainOverNoisedB(:,1+numCPE_tot:end));   
        params.D = D_cell;
        params.R_gNB = R_gNB(:,:,:,1+numCPE_tot:end);
        params.R_cpe = [];
        params.R_interue = [];
        params.R_ue = R_ue; 
        nbrOfRealizations = params.nbrOfRealizations;
        rate_dl = zeros(K,nbrOfRealizations);
        for n = 1:nbrOfRealizations
            [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
            % rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
            rate_dl(:,n) = compute_link_rates_OFDM(params, channel_dl);                                              
        end
        mean_rate_dl_cell = mean(rate_dl,2);
        for idxSI = 1:length(SI_cancel_arr)
            params.SI_cancel_factor = 10^(0.1*SI_cancel_arr(idxSI)); %inter-BS interference cancellation factor
            params.numUE = 0;
            params.numCPE = numCPE_tot;
            K_FWA = params.numCPE;
            K = params.numCPE + M_sectors*params.numUE; 
            params.CPE_locations = CPE_locations;
            params.BETA = db2pow(gainOverNoisedB(:,1:numCPE_tot));   
            params.D = D_FWA;
            params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
            params.R_cpe = R_cpe;
            params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
            params.R_ue = []; 
            rate_dl = zeros(K,nbrOfRealizations);
            for n = 1:nbrOfRealizations
                [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
            end
            mean_rate_dl_FWA = mean(rate_dl,2);
            save_old_mean_FWA = mean_rate_dl_FWA;
            Band_FWA = params.Band;
            % for idxrmin = 1:length(r_min_arr)
            %     sum_FWA_se = 0;
            %     params.r_min_FWA = r_min_arr(idxrmin);
            %     K_FWA_max = 0;
            %     params.Band = Band_FWA;
            %     if (params.Band > 0)
            %         while any(mean_rate_dl_FWA > params.r_max_FWA)
            %             CPE_idxs = find(mean_rate_dl_FWA > params.r_max_FWA);
            %             K_FWA_max = K_FWA_max + numel(CPE_idxs);
            %             params.numCPE = params.numCPE - numel(CPE_idxs);
            %             sum_FWA_se = sum_FWA_se + sum(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA))/params.Band;
            %             params.Band = params.Band*(1-(params.r_min_FWA/min(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA))));
            %             params.R_gNB(:,:,:,CPE_idxs) = [];
            %             params.R_cpe(:,:,:,CPE_idxs) = [];
            %             params.D(:,CPE_idxs) = [];
            %             rate_dl = zeros(K-K_FWA_max,nbrOfRealizations);
            %             for n = 1:nbrOfRealizations
            %                 [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
            %                 rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
            %             end
            %             mean_rate_dl_FWA = mean(rate_dl,2);
            %         end
            %         [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
            %         CPE_idxs = find(FWA_util == 0);
            %         K_FWA_max = K_FWA_max + sum(FWA_util>0);
            %         rate_dl = zeros(sum(FWA_util > 0),nbrOfRealizations);
            %         params.numCPE = params.numCPE - numel(CPE_idxs);
            %         params.R_gNB(:,:,:,CPE_idxs) = [];
            %         params.R_cpe(:,:,:,CPE_idxs) = [];
            %         params.D(:,CPE_idxs) = [];
            %         for n = 1:nbrOfRealizations
            %             [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
            %             rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                          
            %         end
            %         mean_rate_dl_FWA = mean(rate_dl,2);
            %         sum_FWA_se = sum_FWA_se + sum(mean_rate_dl_FWA)/params.Band;
            %         params.Band = params.Band*(1-params.r_min_FWA/min(mean_rate_dl_FWA));
            %     end
            %     params.CPE_locations = CPE_locations;
            %     params.numCPE = numCPE_tot;
            %     params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
            %     params.R_cpe = R_cpe;
            %     params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
            %     params.R_ue = []; 
            %     params.D = D_FWA;
            %     mean_rate_dl_FWA = save_old_mean_FWA;
            %     %% Recording the Results
            % 
            %     %Taking care of folder directory creation etc
            %     dataFolder = 'resultData';
            %     rateFolder = strcat(dataFolder,'/FWA_SE_comp_multi_cell_present_results');
            %     if not(isfolder(dataFolder))
            %         mkdir(dataFolder)
            %     end
            %     if not(isfolder(rateFolder))
            %         mkdir(rateFolder)
            %     end
            % 
            % 
            %     %Saving all results as a structure
            %     dataDescription = {'simOutputs is a 4D array';...
            %         ', for mesh of params ordered as follows';...
            %         'First Dimension: discovery_time';...
            %         'Second Dimension: FailureDetectionTime';...
            %         'Third Dimension: connection_time (RACH)';...
            %         'Fourth Dimension: signalingAfterRachTime';...
            %         '=================================';...
            %         'Each element is a struct'};
            % 
            %     deployRange = params.deployRange;
            %     numBS = size(params.locationsBS,1);
            %     result_string = strcat('/results_numFWA_',num2str(params.SC), 'SC_', num2str(numCPE_all),...
            %         'CPE_',num2str(lambda_BS(idxBSDensity)),...
            %         'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
            %         'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
            %         'lambdaUE_', num2str(deployRange),'deployRange_', ...
            %         num2str(params.r_min_FWA/10^6),'min_FWA_rate', aID);
            %     % result_string = strcat('/results_numFWA_',num2str(params.SC), 'SC_', num2str(numCPE_all),...
            %     %     'CPE_',num2str(lambda_BS(idxBSDensity)),...
            %     %     'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
            %     %     'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
            %     %     'lambdaUE_', num2str(deployRange),'deployRange_', ...
            %     %     num2str(params.r_min_FWA/10^6),'min_FWA_rate', string(aID));
            %     recording_text_file_string = strcat(rateFolder,result_string,'.csv');
            %     fileID = fopen(recording_text_file_string,'w');
            %     output_categories = ['lambdaBS,','lambdaSC,','numCPE,','lambdaUE,',...
            %     'deployRange,','r_min_cell,','r_min_FWA,','max_FWA,', 'Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
            %     fprintf(fileID,output_categories);
            %     formatSpec = '%d,%d,%d,%d,%d,%f,%f,%d,%f,%f,%f,%f\n';
            %     fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_SC(idxBSDensity),numCPE_all, ...
            %     lambda_UE(idxUEDensity),deployRange,params.r_min_cell,params.r_min_FWA,K_FWA_max,Band,Band_FWA, mean(mean_rate_dl_cell)/Band, sum_FWA_se/K_FWA_max);
            %     fclose(fileID);
            % end
            %Taking care of folder directory creation etc
            dataFolder = 'resultData';
            rateFolder = strcat(dataFolder,'/FWA_SE_comp_multi_cell_no_rmin_SI_cancel');
            if not(isfolder(dataFolder))
                mkdir(dataFolder)
            end
            if not(isfolder(rateFolder))
                mkdir(rateFolder)
            end
    
    
            %Saving all results as a structure
            dataDescription = {'simOutputs is a 4D array';...
                ', for mesh of params ordered as follows';...
                'First Dimension: discovery_time';...
                'Second Dimension: FailureDetectionTime';...
                'Third Dimension: connection_time (RACH)';...
                'Fourth Dimension: signalingAfterRachTime';...
                '=================================';...
                'Each element is a struct'};
    
            deployRange = params.deployRange;
            numBS = size(params.locationsBS,1);
            result_string = strcat('/results_numFWA_',num2str(params.SC), ...
                'SC_', num2str(numCPE_all),...
                'CPE_',num2str(lambda_BS(idxBSDensity)),...
                'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
                'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
                'lambdaUE_', num2str(deployRange), ...
                'deployRange_',num2str(params.SI_cancel_factor), ...
                'SI_cancel_',aID);
            % result_string = strcat('/results_numFWA_',num2str(params.SC), 'SC_', num2str(numCPE_all),...
            %     'CPE_',num2str(lambda_BS(idxBSDensity)),...
            %     'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
            %     'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
            %     'lambdaUE_', num2str(deployRange),'deployRange_', string(aID));
            recording_text_file_string = strcat(rateFolder,result_string,'.csv');
            fileID = fopen(recording_text_file_string,'w');
            output_categories = ['lambdaBS,','lambdaSC,','numCPE,','lambdaUE,',...
            'deployRange,','SIcancel,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n'];
            fprintf(fileID,output_categories);
            formatSpec = '%d,%d,%d,%d,%d,%d,%f,%f,%f,%f\n';
            fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_SC(idxBSDensity),numCPE_all, ...
            lambda_UE(idxUEDensity),deployRange,SI_cancel_arr(idxSI),Band,Band_FWA, mean(mean_rate_dl_cell)/Band, mean(mean_rate_dl_FWA)/Band_FWA);
            fclose(fileID);
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd) 
% end