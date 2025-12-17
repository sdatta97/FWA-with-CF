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
    aID = '0099';
end
% for aID = 1:99
%RNG seed.
rng(str2double(aID),'twister');
% rng(aID,'twister');
%% GUE channel parameters
params.K_Factor = 9;         %dB -- %rician factor Ground UE  % if beta_gains=1
params.RAYLEIGH=0;   %1= rayleigh, % 0=rician
params.Perf_CSI =1;
params.cov_area = 1; %0.25; % 4; %km
%% Mobility params
params.MOBILE = 1;
params.ue_velocity = 11.176; %25 mph
params.Ts = 1.67e-5; %10us
params.fc = 3.7e9; %7.9 GHz FR3
params.c = 3e8; %speed of light
params.coherence_time = 1e-3;
params.n_samples = floor(params.coherence_time/params.Ts);
% params.mob_rho = (besselj(0,2*pi*params.Ts*params.ue_velocity*params.fc/params.c))^params.n_samples;
params.mob_rho = besselj(0,2*pi*params.Ts*params.ue_velocity*params.fc/params.c);
%% Economics params
params.price_cell_per_bit = 4;
params.price_FWA_per_bit = 0.14;
params.r_min_cell = 35e6;
params.r_max_FWA = 100e6;
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

numCPE_all = 50; %5:5:20;
Band = 100e6;
%Prepare to save simulation results
%% Room Setup, UE placement, UE height
params.deployRange = 200; %20:20:100;
params.coverageRange_sub6 = 430;
params.num_antennas_per_gNB = 64;
params.num_antennas_per_sc = 16;
params.rho_tot = 10^(0.1*75)*(Band/1e8); %10^(3.6)*params.num_antennas_per_gNB; %200;
params.rho_tot_sc = 10^(0.1*55); %10^(3.6)*params.num_antennas_per_sc;
%Number of antennas per UE
params.N_UE_FWA = 8;
params.N_UE_cell = 1; %4;
params.hr = 1;
params.hr_cpe = 3;
params.ht_bs = 15;
params.ht_sc = 5;
lambda_BS = 5; %([5 6 7 8 9 10]).^2;
lambda_SC = 0; %([5 6 7 8 9 10]).^2;
lambda_UE = 1000; %0:250:1000; %200:10:250; %150; %100:50:200; %[30:20:90, 100]; %100;
params.Lmax = 3;
params.preLogFactor = 1;
params.loss_pc_FWA = 5;
%Number of channel realizations per setup
params.nbrOfRealizations = 100;

%% UE angular coverage range (full 360 coverage for now)
lookAngleCell{1} = [0,360];
r_min_arr = 90e6; %1e6*(10:20:90);
%% Simulation FR1 setup
%% CPE locations
RCPE =  params.deployRange*sqrt(rand(numCPE_all,1)); %location of UEs (distance from origin)
angleCPE = 2*pi*rand(numCPE_all,1);%location of UEs (angle from x-axis)
CPE_locations = [RCPE.*cos(angleCPE), RCPE.*sin(angleCPE)];
for idxBSDensity = 1:length(lambda_BS)
    %% gNB locations
    params.numGNB = ceil(lambda_BS(idxBSDensity)*pi*(params.deployRange/1000)^2);
    params.RgNB = 0; %params.deployRange*sqrt(rand(params.numGNB,1));
    params.angleGNB = 2*pi*rand(params.numGNB,1);
    params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
    for idxUEDensity = 1:length(lambda_UE)
        %% UE locations
        params.numUE = ceil(lambda_UE(idxUEDensity)*pi*(params.deployRange/1000)^2);
        RUE = params.deployRange*sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
        angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
        params.UE_locations = [RUE.*cos(angleUE), RUE.*sin(angleUE)];  
        params.numCPE = numCPE_all;
        params.CPE_locations = CPE_locations;
        params.Band = Band; %Communication bandwidth
        [gainOverNoisedB,R_gNB,R_cpe,R_ue,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,str2double(aID));
        % [gainOverNoisedB,R_gNB,R_cpe,R_ue,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,aID);
        K_FWA_max = 0;
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
        K = params.numCPE + params.numUE; 
        params.BETA = db2pow(gainOverNoisedB(:,1+numCPE_all:end));   
        if params.SC
            params.D = D_small(:,1+numCPE_all:end);
        else
            params.D = D(:,1+numCPE_all:end);
        end
        params.R_gNB = R_gNB(:,:,:,1+numCPE_all:end);
        params.R_cpe = [];
        params.R_ue = R_ue; 
        nbrOfRealizations = params.nbrOfRealizations;
        rate_dl = zeros(K,nbrOfRealizations);
        for n = 1:nbrOfRealizations
            [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
            rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
        end
        mean_rate_dl_cell = mean(rate_dl,2);
        UE_idxs = find(mean_rate_dl_cell./params.Band > 2*(params.r_min_cell./params.Band));
        if (~isempty(UE_idxs))
            UE_locations = params.UE_locations;
            numUE = params.numUE;
            params.UE_locations = params.UE_locations(UE_idxs,:);
            params.numUE = size(params.UE_locations,1);
            K_FWA = params.numCPE;
            K = params.numCPE + params.numUE; 
            params.BETA = db2pow(gainOverNoisedB(:,numCPE_all+UE_idxs));   
            if params.SC
                params.D = D_small(:,numCPE_all+UE_idxs);
            else
                params.D = D(:,numCPE_all+UE_idxs);
            end
            params.R_gNB = R_gNB(:,:,:,numCPE_all+UE_idxs);
            params.R_cpe = [];
            params.R_ue = R_ue(:,:,:,UE_idxs); 
            nbrOfRealizations = params.nbrOfRealizations;
            rate_dl = zeros(K,nbrOfRealizations);
            for n = 1:nbrOfRealizations
                [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
            end
            mean_rate_dl_cell_new = mean(rate_dl,2);
            band_allotted = params.r_min_cell/min(mean_rate_dl_cell_new./params.Band);
            params.Band = params.Band - band_allotted;
            params.UE_locations = UE_locations;
            params.UE_locations(UE_idxs,:) = [];                
            params.numUE = numUE - numel(UE_idxs);
            if (params.numUE > 0)
                K_FWA = params.numCPE;
                K = params.numCPE + params.numUE; 
                UE_old_idxs = UE_idxs;
                UE_idxs = 1:1:numUE;
                UE_idxs(UE_old_idxs) = [];
                params.BETA = db2pow(gainOverNoisedB(:,numCPE_all+UE_idxs));   
                if params.SC
                    params.D = D_small(:,numCPE_all+UE_idxs);
                else
                    params.D = D(:,numCPE_all+UE_idxs);
                end
                params.R_gNB = R_gNB(:,:,:,numCPE_all+UE_idxs);
                params.R_cpe = [];
                params.R_ue = R_ue(:,:,:,UE_idxs); 
                nbrOfRealizations = params.nbrOfRealizations;
                rate_dl = zeros(K,nbrOfRealizations);
                for n = 1:nbrOfRealizations
                    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                    rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                end
                mean_rate_dl_cell = mean(rate_dl,2);
                if all(mean_rate_dl_cell>=params.r_min_cell)
                    params.Band = max(params.Band*(1-params.r_min_cell/min(mean_rate_dl_cell)),0);
                    if (params.Band > 0)
                        mean_rate_dl_cell = mean_rate_dl_cell*params.r_min_cell/min(mean_rate_dl_cell); 
                    end
                else
                    params.Band = 0;
                end
            end
        end
        params.numUE = 0;
        params.numCPE = numCPE_all;
        K_FWA = params.numCPE;
        K = params.numCPE + params.numUE; 
        params.CPE_locations = CPE_locations;
        params.BETA = db2pow(gainOverNoisedB(:,1:numCPE_all));   
        if params.SC
            params.D = D_small(:,1:numCPE_all);
        else
            params.D = D(:,1:numCPE_all);
        end
        params.R_gNB = R_gNB(:,:,:,1:numCPE_all);
        params.R_cpe = R_cpe;
        params.R_ue = []; 
        rate_dl = zeros(K,nbrOfRealizations);
        for n = 1:nbrOfRealizations
            [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
            rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
        end
        mean_rate_dl_FWA = mean(rate_dl,2);
        for idxrmin = 1:length(r_min_arr)
            params.r_min_FWA = r_min_arr(idxrmin);
            K_FWA_max = 0;
            if (params.Band > 0)
                new_CPE_idxs = find(mean_rate_dl_FWA./params.Band > 2*(params.r_min_FWA./params.Band));
                K_FWA_max = K_FWA_max + numel(new_CPE_idxs);
                if (~isempty(new_CPE_idxs))
                    CPE_locations = params.CPE_locations;
                    numCPE = params.numCPE;
                    params.CPE_locations = params.CPE_locations(new_CPE_idxs,:);
                    params.numCPE = size(params.CPE_locations,1);
                    K_FWA = params.numCPE;
                    K = params.numCPE + params.numUE; 
                    params.BETA = db2pow(gainOverNoisedB(:,new_CPE_idxs));   
                    if params.SC
                        params.D = D_small(:,new_CPE_idxs);
                    else
                        params.D = D(:,new_CPE_idxs);
                    end
                    params.R_gNB = R_gNB(:,:,:,new_CPE_idxs);
                    params.R_cpe = R_cpe(:,:,:,new_CPE_idxs);
                    params.R_ue = []; 
                    nbrOfRealizations = params.nbrOfRealizations;
                    rate_dl = zeros(K,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                    end
                    mean_rate_dl_FWA_new = mean(rate_dl,2);
                    band_allotted = params.r_min_FWA/min(mean_rate_dl_FWA_new./params.Band);
                    params.Band = params.Band - band_allotted;
                end
                params.CPE_locations = CPE_locations;
                params.CPE_locations(new_CPE_idxs,:) = [];                
                params.numCPE = numCPE - numel(new_CPE_idxs);
                idxs_new = [];
                if (params.numCPE > 0)
                    K_FWA = params.numCPE;
                    K = params.numCPE + params.numUE; 
                    CPE_old_idxs = new_CPE_idxs;
                    CPE_idxs = 1:1:numCPE;
                    CPE_idxs(CPE_old_idxs) = [];
                    params.BETA = db2pow(gainOverNoisedB(:,CPE_idxs));   
                    if params.SC
                        params.D = D_small(:,CPE_idxs);
                    else
                        params.D = D(:,CPE_idxs);
                    end
                    params.R_gNB = R_gNB(:,:,:,CPE_idxs);
                    params.R_cpe = R_cpe(:,:,:,CPE_idxs);
                    params.R_ue = []; 
                    nbrOfRealizations = params.nbrOfRealizations;
                    rate_dl = zeros(K,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                    idxs_new = find(FWA_util>0);
                    K_FWA_max = K_FWA_max + numel(idxs_new);
                    params.Band = max(params.Band - params.r_min_FWA/(min(mean_rate_dl_FWA)./params.Band),0);                    
                    params.CPE_locations(idxs_new,:) = [];                
                    params.numCPE = params.numCPE - numel(idxs_new);
                end
                if ((params.Band > 0) && (params.numCPE>0))
                    K_FWA = params.numCPE;
                    K = params.numCPE + params.numUE; 
                    CPE_idxs = 1:1:numCPE;
                    CPE_idxs(CPE_old_idxs) = [];
                    CPE_idxs(idxs_new) = [];
                    params.BETA = db2pow(gainOverNoisedB(:,CPE_idxs));   
                    if params.SC
                        params.D = D_small(:,CPE_idxs);
                    else
                        params.D = D(:,CPE_idxs);
                    end
                    params.R_gNB = R_gNB(:,:,:,CPE_idxs);
                    params.R_cpe = R_cpe(:,:,:,CPE_idxs);
                    params.R_ue = []; 
                    nbrOfRealizations = params.nbrOfRealizations;
                    rate_dl = zeros(K,nbrOfRealizations);
                    for n = 1:nbrOfRealizations
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                    K_FWA_max = K_FWA_max + numel(find(FWA_util>0));
                end
            end
            params.CPE_locations = CPE_locations;
            params.numCPE = numCPE;
            params.R_gNB = R_gNB(:,:,:,1:numCPE_all);
            params.R_cpe = R_cpe;
            params.R_ue = []; 
            %% Recording the Results
        
            %Taking care of folder directory creation etc
            dataFolder = 'resultData';
            rateFolder = strcat(dataFolder,'/FWA_cell_best_results');
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
            result_string = strcat('/results_numFWA_',num2str(params.SC), 'SC_', num2str(numCPE_all),...
                'CPE_',num2str(lambda_BS(idxBSDensity)),...
                'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
                'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
                'lambdaUE_', num2str(deployRange),'deployRange_', ...
                num2str(params.r_min_FWA/10^6),'min_FWA_rate', aID);
            % result_string = strcat('/results_numFWA_',num2str(params.SC), 'SC_', num2str(numCPE_all),...
            %     'CPE_',num2str(lambda_BS(idxBSDensity)),...
            %     'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
            %     'lambdaSC_',num2str(lambda_UE(idxUEDensity)),...
            %     'lambdaUE_', num2str(deployRange),'deployRange_', ...
            %     num2str(params.r_min_FWA/10^6),'min_FWA_rate', string(aID));
            recording_text_file_string = strcat(rateFolder,result_string,'.csv');
            fileID = fopen(recording_text_file_string,'w');
            output_categories = ['lambdaBS,','lambdaSC,','numCPE,','lambdaUE,',...
            'deployRange,','r_min_FWA,','max_FWA\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
            fprintf(fileID,output_categories);
            formatSpec = '%d,%d,%d,%d,%d,%f,%d,%d\n';
            fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_SC(idxBSDensity),numCPE_all, ...
            lambda_UE(idxUEDensity),deployRange,params.r_min_FWA,K_FWA_max);
            fclose(fileID);
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd) 
% end