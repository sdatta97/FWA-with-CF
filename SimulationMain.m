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
    aID = '0022';
end
%RNG seed.
rng(str2double(aID),'twister');
% for aID = 1:10
%     RNG seed.
%     rng(aID,'twister');
%% GUE channel parameters
params.K_Factor = 9;         %dB -- %rician factor Ground UE  % if beta_gains=1
params.RAYLEIGH=0;   %1= rayleigh, % 0=rician
params.Perf_CSI =1;
params.cov_area = 1; %0.25; % 4; %km
%%
params.TAU_P_K_by_two = 0; %1;  
params.CH_estimation = 0;  % 1= have channel estimation
%%
params.LB=1;  %Lower bound
params.UB =1;  %Upper bound
params.no_of_rea = 1;     % no.of channel realizations
%%
% snr_db = -50:10:40;
params.snr_db = 40;
params.ASD_VALUE = 0;%[0,0.25,0.5,0.75,1];  % [0,30,10]; %
params.ASD_CORR = 0;
params.Kt_Kr_vsUE  = 0; %0.175^2; %0.175^2; %[1,2,3,4];  %to save 1=AP 0.1,UE=0.1;  2=AP 0.1,UE=0.3;  3=AP 0.3,UE=0.1

params.pilot_pow = 100;  % 0.1W   % UL pilot. power (W)
params.noiseFigure = 9; % gue
params.sigma_sf =4;
params.Band = 100e6; %Communication bandwidth

params.SC = 0; %0 - CF, 1 - SC
%% Define simulation setup

%Angular standard deviation in the local scattering model (in radians)
params.ASD_varphi = deg2rad(30); %azimuth angle
params.ASD_theta = 0; %deg2rad(15);  %elevation angle

%Total uplink transmit power per UE (mW)
params.p = 100;

%Total downlink transmit power per AP (mW)
% rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];

%Power factor division
% p_fac_arr = 10; %10.^(0:1:2);
params.p_fac = 10;
numCPE_arr = 10:10:50; %5:5:20;

%Prepare to save simulation results

% rng(2,'twister');
%%
% load('params.mat')
params.simTime = 60*60; %sec Total Simulation time should be more than 100.
%% Room Setup, UE placement, UE height
params.deployRange = 200; %20:20:100;
params.deployRange_sub6 = 1000;
params.coverageRange_sub6 = 430;
params.num_antennas_per_gNB = 64;
params.num_antennas_per_sc = 16;
params.rho_tot = 10^(3.6)*params.num_antennas_per_gNB; %200;
params.rho_tot_sc = 10^(3.6)*params.num_antennas_per_sc;
%Number of antennas per UE
params.N_UE_FWA = 8;
params.N_UE_cell = 1; %4;
params.hr = 1;
params.ht_bs = 15;
params.ht_sc = 5;
lambda_BS = 20; %([5 6 7 8 9 10]).^2;
lambda_UE = 100; %250:250:1000; %200:10:250; %150; %100:50:200; %[30:20:90, 100]; %100;
params.Lmax = 4;
lb_thresh = 0; %0.1; %[0:0.05:0.1 0.5 1];
params.preLogFactor = 1;
%Number of channel realizations per setup
params.nbrOfRealizations = 100;

%% UE angular coverage range (full 360 coverage for now)
lookAngleCell{1} = [0,360];
for idxnumCPE = 1:length(numCPE_arr)
    params.numCPE = numCPE_arr(idxnumCPE);
    %%UE location
    deployRange = params.deployRange; %(idxdeployRange);
    deployRange_sub6 = params.deployRange_sub6; %(idxdeployRange);
    params.RUE =  deployRange*sqrt(rand(params.numCPE,1)); %location of UEs (distance from origin)
    params.angleUE = 2*pi*rand(params.numCPE,1);%location of UEs (angle from x-axis)
    params.CPE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
    %% Simulation FR1 setup
    for idxUEDensity = 1:length(lambda_UE)
        params.numUE = ceil(lambda_UE(idxUEDensity)*pi*(params.deployRange_sub6/1000)^2);
        params.RUE_sub6 = params.deployRange_sub6*sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
        params.angleUE_sub6 = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
        params.UE_locations =  [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
        for idxBSDensity = 1:length(lambda_BS)
            %% gNB locations
            params.numGNB = ceil(lambda_BS(idxBSDensity)*pi*(params.deployRange/1000)^2);
            params.RgNB = params.deployRange*sqrt(rand(params.numGNB,1));
            params.angleGNB = 2*pi*rand(params.numGNB,1);
            params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
            params.numSC = ceil(lambda_BS(idxBSDensity)*pi*((params.deployRange_sub6/1000)^2-(params.deployRange/1000)^2));
            params.Rsc = params.deployRange_sub6*sqrt(rand(params.numSC,1));
            params.anglesc = 2*pi*rand(params.numSC,1);
            params.locationsSC = [params.Rsc.*cos(params.anglesc), params.Rsc.*sin(params.anglesc)];
            
            K = params.numCPE + params.numUE; 
            snr_db = params.snr_db;
            ASD_VALUE = params.ASD_VALUE;
            ASD_CORR = params.ASD_CORR;
            Kt_Kr_vsUE = params.Kt_Kr_vsUE;
            K_Factor = params.K_Factor;
            RAYLEIGH=params.RAYLEIGH;   %1= rayleigh, % 0=rician
            Perf_CSI = params.Perf_CSI;
            cov_area = params.cov_area;
            %%
            no_of_rea = params.no_of_rea;     % no.of channel realizations
            %%
            pilot_pow = params.pilot_pow; 
            noiseFigure = params.noiseFigure;
            sigma_sf = params.sigma_sf;
            Band = params.Band; %Communication bandwidth
            [gainOverNoisedB,R_gNB,R_cpe,R_ue,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,str2double(aID));
            [gainOverNoisedB_sc,R_sc,R_cpe_sc,R_ue_sc,D_sc,D_small_sc,SCpositions,distances_sc] = generateSetup_SC(params,str2double(aID),pilotIndex,D,D_small,UEpositions);
            params.BETA = db2pow(gainOverNoisedB);   
            params.BETA_sc = db2pow(gainOverNoisedB_sc);   
            % params.D = D;
            if params.SC
                params.D = D_small;
                params.D_sc = D_small_sc;
            else
                params.D = D;
                params.D_sc = D_sc;
            end
            params.R_gNB = R_gNB;
            params.R_sc = R_sc;
            params.R_cpe = R_cpe;
            params.R_ue = R_ue; 
            params.R_cpe_sc = R_cpe_sc;
            params.R_ue_sc = R_ue_sc; 
            nbrOfRealizations = params.nbrOfRealizations;
            rate_dl = zeros(K,nbrOfRealizations);
            for n = 1:nbrOfRealizations
                [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                [channel_dl_sc, channel_est_dl_sc,channel_dl_FWA_sc, channel_est_dl_FWA_sc] = computePhysicalChannels_sub6_MIMO_sc(params);
                rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_dl_sc, channel_est_dl, channel_est_dl_sc, channel_dl_FWA, channel_dl_FWA_sc, channel_est_dl_FWA, channel_est_dl_FWA_sc);                                              
            end
            mean_rate_dl = mean(rate_dl,2);
            %% Recording the Results
        
            %Taking care of folder directory creation etc
            dataFolder = 'resultData';
            rateFolder = strcat(dataFolder,'/FWA_cell_revised_results');
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
    
            numCPE = params.numCPE;
            numUE = params.numUE;
            numBS = size(params.locationsBS,1);
            result_string = strcat('/results_',num2str(params.SC), 'SC_', num2str(numCPE_arr(idxnumCPE)),...
                'CPE_',num2str(lambda_BS(idxBSDensity)),...
                'lambdaBS_',num2str(lambda_UE(idxUEDensity)),...
                'lambdaUE_',num2str(deployRange),...
            'deployRange_',aID);
            recording_text_file_string = strcat(rateFolder,result_string,'.csv');
            fileID = fopen(recording_text_file_string,'w');
            output_categories = ['lambdaBS,','numCPE,','lambdaUE,',...
                'deployRange,', 'mean_rate_FWA,','mean_rate_cell\n'];
            fprintf(fileID,output_categories);

            mean_rate_FWA = mean(mean_rate_dl(1:numCPE));
            mean_rate_cell = mean(mean_rate_dl((1+numCPE):end));
            
            formatSpec = '%d,%d,%d,%f,%f,%f,%f,%.16f,%.16f,%.16f,%.16f,%.16f\n';
            fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),numCPE,lambda_UE(idxUEDensity),...
            deployRange,mean_rate_FWA,mean_rate_cell);
            fclose(fileID);
        end
    end
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)
% end
