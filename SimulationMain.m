close all;
clear;
tStart = tic;
% aID = getenv('SLURM_ARRAY_TASK_ID');

% This is for running on a cluster in parallel
% the bash script should give the aID as input
% if (isempty(aID))
%     warning('aID is empty. Trying SLURM ID.')
%     aID = getenv('SLURM_ARRAY_TASK_ID');
% end
% if(isempty(aID))
%     warning('aID is empty. Replacing it with 0010.')
%     aID = '0022';
% end
% %RNG seed.
% rng(str2double(aID),'twister');
for aID = 0:99
    %RNG seed.
    rng(aID,'twister');
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
    
    
    %% Define simulation setup
    
    %Angular standard deviation in the local scattering model (in radians)
    params.ASD_varphi = deg2rad(30); %azimuth angle
    params.ASD_theta = 0; %deg2rad(15);  %elevation angle
    
    %Total uplink transmit power per UE (mW)
    params.p = 100;
    
    %Total downlink transmit power per AP (mW)
    % rho_tot_arr = [10:10:100, 200:100:1000, 2000:1000:10000];
    
    %Power factor division
    p_fac_arr = 10; %10.^(0:1:2);
    % params.p_fac = 10;
    percent_fr2_UE_arr = 5; %5:5:20;
    
    %Prepare to save simulation results
    
    % rng(2,'twister');
    %%
    % load('params.mat')
    params.simTime = 60*60; %sec Total Simulation time should be more than 100.
    %% Room Setup, UE placement, UE height
    % We are considering an outdoor scenario where the UE is located at the
    % center and gNBs are distributed around the UE. We only need to consider
    % the coverageRange amount of distance from the UE.
    params.deployRange = 300; %20:20:100;
    params.deployRange_sub6 = 1000;
    
    params.coverageRange_sub6 = 430;
    length_area_sub6 = 2*params.coverageRange_sub6;   
    width_area_sub6 = 2*params.coverageRange_sub6;
    % length_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);   
    % width_area_sub6 = 2*(params.deployRange + params.coverageRange_sub6);
    height_transmitter = 5;
    height_transmitter_sub6 = 4; % 5;
    params.areaDimensions_sub6 = [width_area_sub6, length_area_sub6, height_transmitter_sub6];
    
    params.hr = 1.4; %height receiver (UE), approximately the height a human holds the phone
    params.ht = height_transmitter; %height transmitter (BS)
    params.ht_sub6 = height_transmitter_sub6; %height transmitter (BS)
    params.num_antennas_per_gNB = 64;
    params.rho_tot = 10^(3.6)*params.num_antennas_per_gNB; %200;
    
    % params.num_antennas_per_gNB = 8;
    %Number of antennas per UE
    params.N_UE_FWA = 8;
    params.N_UE_cell = 1; %4;
    rmin_arr = 4*10^8;
    lambda_BS = 5; %([5 6 7 8 9 10]).^2;
    lambda_UE_sub6 = 250; %250:250:1000; %200:10:250; %150; %100:50:200; %[30:20:90, 100]; %100;
    params.Lmax = 2;
    % for idxnumUEsub6 = 1:length(numUE_sub6_arr)
    lb_thresh = 0.1; %[0:0.05:0.1 0.5 1];
    for idxnumUE = 1:length(percent_fr2_UE_arr)
        % n = poissrnd((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
        % while(n==0)
        %     n = poissrnd((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
        % end
        % params.numUE = n;
        % params.numUE = ceil((percent_fr2_UE_arr(idxnumUE)/100)*lambda_UE_sub6*pi*(params.deployRange/1000)^2);
        params.numUE = 20;
        %%UE location
        deployRange = params.deployRange; %(idxdeployRange);
        deployRange_sub6 = params.deployRange_sub6; %(idxdeployRange);
        params.RUE =  deployRange_sub6*sqrt(rand(params.numUE,1)); %location of UEs (distance from origin)
        params.angleUE = 2*pi*rand(params.numUE,1);%location of UEs (angle from x-axis)
        params.UE_locations = [params.RUE.*cos(params.angleUE), params.RUE.*sin(params.angleUE)];
        %% Simulation FR1 setup
        for idxUEDensity = 1:length(lambda_UE_sub6)
            params.ue_rearranged = [];        
            %%UE locations
            % n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
            % while (n==0)
            %     n = poissrnd(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);       
            % end
            n = ceil(lambda_UE_sub6(idxUEDensity)*pi*(params.coverageRange_sub6/1000)^2);
            params.numUE_sub6 = n;
            % params.RUE_sub6 = params.coverageRange_sub6*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            % params.RUE_sub6 = (norm(params.UE_locations(:,1) - params.UE_locations(:,2))/2+params.coverageRange_sub6)*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            params.RUE_sub6 = (max(sqrt(sum(params.UE_locations.^2,2)))+params.coverageRange_sub6)*sqrt(rand(params.numUE_sub6,1)); %location of UEs (distance from origin)
            params.angleUE_sub6 = 2*pi*rand(params.numUE_sub6,1);%location of UEs (angle from x-axis)
            % params.UE_locations_sub6 =  mean(params.UE_locations,1) + [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
            params.UE_locations_sub6 =  [params.RUE_sub6.*cos(params.angleUE_sub6), params.RUE_sub6.*sin(params.angleUE_sub6)];   
            for idxBSDensity = 1:length(lambda_BS)
                %% gNB locations
                params.numGNB = ceil(lambda_BS(idxBSDensity)*(params.deployRange_sub6/1000)^2);
                % params.xGNB = (-params.deployRange_sub6/2):(params.deployRange_sub6/(sqrt(params.numGNB)-1)):(params.deployRange_sub6/2);
                % params.yGNB = (-params.deployRange_sub6/2):(params.deployRange_sub6/(sqrt(params.numGNB)-1)):(params.deployRange_sub6/2);
                params.RgNB = params.deployRange_sub6*sqrt(rand(params.numGNB,1));
                params.angleGNB = 2*pi*rand(params.numGNB,1);
                params.locationsBS = [params.RgNB.*cos(params.angleGNB), params.RgNB.*sin(params.angleGNB)];
                % params.locationsBS = (combvec(params.xGNB,params.yGNB)).';
                params.coverageRange = 100; %(params.deployRange_sub6/(sqrt(params.numGNB)-1))/sqrt(2);%125*sqrt(2); %100;
                length_area = 2*params.coverageRange;   
                width_area = 2*params.coverageRange;
                % length_area = 2*(params.deployRange + params.coverageRange);   
                % width_area = 2*(params.deployRange + params.coverageRange);
                params.areaDimensions = [width_area, length_area, height_transmitter];
                        %Length of the coherence block
        %                 params.tau_c = 200;
                    
                    %Compute number of pilots per coherence block
        %                 params.tau_p = params.numUE+params.numUE_sub6;
                    
                    %Compute the prelog factor assuming only downlink data transmission
            %         params.preLogFactor = (params.tau_c-params.tau_p)/params.tau_c;
                params.preLogFactor = 1;
        
                %Number of setups with random UE locations
                params.nbrOfSetups = 100;
                        
                      
                %Number of channel realizations per setup
                params.nbrOfRealizations = 100;
                
                %% PHY layer params
                params.scs_mmw = 2e9;     %not using this parameter now
                params.scs_sub6 = [0.5*(params.Band), 0.5*(params.Band)];   %sub-6 GHz bandwidth 100 MHz
                params.num_sc_mmw = 1;    %not using this parameter now
                params.num_sc_sub6 = 2;   %sub-6 GHz considered as one full band                    
                %% UE angular coverage range (full 360 coverage for now)
                lookAngleCell{1} = [0,360];
                
                %% Blocker Properties and Simulation Duration
                params.lambdaBlockers = 0.01; %How many blockers around
                % params.numBlockers = 4*(params.coverageRange)^2*params.lambdaBlockers;
                params.numBlockers = floor(pi*(params.coverageRange)^2*params.lambdaBlockers);
        %         params.numBlockers = 4*(params.coverageRange_sub6)^2*params.lambdaBlockers;
                params.V = 1; %velocity of blocker m/s
                % 160-190 cm truncated gaussian with mean at 3 sigma to each sides.
                % params.hb = (175 + TruncatedGaussian(5, [-15,15], [1 params.numBlockers])) / 100;
                params.hb = 1.8*ones(1,params.numBlockers); %height blocker
                params.mu = 2; %Expected bloc dur =1/mu sec
                numUE = params.numUE;
                numUE_sub6 = params.numUE_sub6;
                
                K = params.numUE + params.numUE_sub6;  % --Ground UEs
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
                num_sc_sub6 = params.num_sc_sub6;
                % params.user_sc_alloc = randi([1,num_sc_sub6],K,1);
                params.user_sc_alloc = zeros(K,num_sc_sub6);
    %                 tau_c = params.tau_c;      % coherence block length 
                % [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,str2double(aID));
                [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,1,aID);
                num_sc_sub6 = params.num_sc_sub6;
                params.BETA = db2pow(gainOverNoisedB);   
                params.D = D;
                % params.D = D_small;
                params.R_gNB = R_gNB;
                params.R_ue_mmW = R_ue_mmW;
                params.R_ue_sub6 = R_ue_sub6;                       
                for idxrmin = 1:length(rmin_arr)
                    for idx_p = 1:length(p_fac_arr)
                        for idxlbthres = 1:length(lb_thresh)
                            params.ue_rearranged = [];
                            lb_thres = lb_thresh(idxlbthres);
                            rmin = rmin_arr(idxrmin);
                            rmin_sub6 = 35e6;
                            params.loss_pc_thresh = lb_thres;
                            params.r_min_sub6 = rmin_sub6;  %stores min rate requirement for all sub-6 users
                            params.r_min = rmin;  %stores min rate requirement for all mmWave users
                            params.rate_reduce_threshold = 5e7;
                            params.p_fac = p_fac_arr(idx_p);
                            params.p_fac_rearrange = 1; % 0.1*p_fac_arr(idx_p);  
                            sub6ConnectionState = zeros(params.numUE,1);
                            ue_idx = 1:10;
                            sub6ConnectionState(ue_idx) = 1;
                            [channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW] = computePhysicalChannels_sub6_MIMO(params);
                            [~, ue_idxs_affected] = AP_reassign(params,ue_idx);
                            ue_rearranged = union(ue_idxs_affected, params.ue_rearranged);
                            rate_dl_before_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,zeros(params.numUE,1));                                              
                            lb = quantile(rate_dl_before_handoff(ue_rearranged)./params.Band,params.loss_pc_thresh);
                            bw_alloc = Band - rmin_sub6/lb;
                            if (bw_alloc < 0) %|| isnan(bw_alloc)
                                bw_alloc = 0;
                                params.p_fac = 1;
                                rate_dl_after_handoff = rate_dl_before_handoff;
                            elseif isnan(bw_alloc)
                                rate_dl_after_handoff = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);                                              
                                lb = quantile(rate_dl_after_handoff((1+numUE):(numUE+numUE_sub6)),params.loss_pc_thresh);                                
                            else 
                                params.ue_rearranged = ue_rearranged;
                                ues_not_affected = setdiff((1+numUE):(numUE+numUE_sub6),params.ue_rearranged);
                                params.scs_sub6(1) = bw_alloc;
                                params.scs_sub6(2) = Band - bw_alloc;
                                % user_sc_alloc = ones(numUE+numUE_sub6,params.num_sc_sub6);                               
                                user_sc_alloc = params.user_sc_alloc; %zeros(numUE+numUE_sub6,1);     
                                user_sc_alloc(find(sub6ConnectionState),1) = 1;
                                user_sc_alloc(find(sub6ConnectionState),2) = 0;
                                user_sc_alloc(ues_not_affected,1) = 1;
                                user_sc_alloc(ues_not_affected,2) = 1;
                                user_sc_alloc(params.ue_rearranged,1) = 0;
                                user_sc_alloc(params.ue_rearranged,2) = 1;
                                params.user_sc_alloc = user_sc_alloc;
                                ues_sharing = union(((1:numUE).*sub6ConnectionState),ues_not_affected);
        %                         rate_dl_after_handoff = compute_link_rates_MIMO(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);                                              
                                rate_dl_after_handoff = compute_link_rates_MIMOv4(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState);    
                                ues_not_affected = setdiff((1+numUE):(numUE+numUE_sub6),params.ue_rearranged);
                                lb = quantile(rate_dl_after_handoff(ues_not_affected),params.loss_pc_thresh);
                            end             
                            %% Recording the Results
                        
                            %Taking care of folder directory creation etc
                            dataFolder = 'resultData';
                            rateFolder = strcat(dataFolder,'/ratecomparisonResults');
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
                    
                            numUE = params.numUE;
                            numUE_sub6 = params.numUE_sub6;
                            numBS = size(params.locationsBS,1);
                            min_rate_req = params.r_min;
                            p_fac = params.p_fac;
                            result_string = strcat('/results_',num2str(percent_fr2_UE_arr(idxnumUE)),...
                                'percentfr2UE_',num2str(lambda_BS(idxBSDensity)),...
                                'lambdaBS_',num2str(lambda_UE_sub6(idxUEDensity)),...
                                'lambdaUE_',num2str(deployRange),...
                            'deployRange_',num2str(aID),'Min_rate', num2str(rmin), "Pow_fac", num2str(p_fac), "lb_thres", num2str(100*lb_thres));
                            recording_text_file_string = strcat(rateFolder,result_string,'.csv');
                            fileID = fopen(recording_text_file_string,'w');
                            output_categories = ['UE idx,','lambdaBS,','lambdaUE,',...
                                'deployRange,','minRatereq,','powerFac,','lower_bound_thresh,', 'mean_rate_before_handoff_affected,','mean_rate_after_handoff_affected,','mean_rate_before_handoff_not_affected,','mean_rate_after_handoff_not_affected,','rate_after_handoff_fr2\n'];
                            fprintf(fileID,output_categories);
    
                            mean_rate_before_handoff_affected = mean(rate_dl_before_handoff(params.ue_rearranged));
                            mean_rate_after_handoff_affected = mean(rate_dl_after_handoff(params.ue_rearranged));
                            mean_rate_before_handoff_not_affected = mean(rate_dl_before_handoff(ues_not_affected));
                            mean_rate_after_handoff_not_affected = mean(rate_dl_after_handoff(ues_not_affected));
                            fr2_rate = rate_dl_after_handoff(ue_idx);
                            
                            formatSpec = '%d,%d,%d,%f,%f,%f,%f,%.16f,%.16f,%.16f,%.16f,%.16f\n';
                            fprintf(fileID,formatSpec,ue_idx, lambda_BS(idxBSDensity),lambda_UE_sub6(idxUEDensity),...
                            deployRange, min_rate_req, p_fac, lb_thres, ...
                            mean_rate_before_handoff_affected,mean_rate_after_handoff_affected,mean_rate_before_handoff_not_affected,mean_rate_after_handoff_not_affected,fr2_rate);
                            fclose(fileID);
                    
                        end
                    end
                end
            end
        end
    end
    tEnd = toc(tStart);
    fprintf('Total runtime: %f seconds\n',tEnd)
end
