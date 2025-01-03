function [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,pilotIndex,D,D_small,APpositions,UEpositions,distances] = generateSetup(params,nbrOfSetups,seed)
%This function generates realizations of the simulation setup described in
%Section 5.3.
%
%INPUT:
%L               = Number of APs per setup
%K               = Number of UEs in the network
%N               = Number of antennas per AP
%tau_p           = Number of orthogonal pilots
%nbrOfSetups     = Number of setups with random UE and AP locations
%seed            = Seed number of pseudorandom number generator
%ASD_varphi      = Angular standard deviation in the local scattering model 
%                  for the azimuth angle (in radians)
%ASD_theta       = Angular standard deviation in the local scattering model
%                  for the elevation angle (in radians)
%
%OUTPUT:
%gainOverNoisedB = Matrix with dimension L x K x nbrOfSetups where
%                  element (l,k,n) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k,n) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k,n)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for small-cell setup
%APpositions     = Vector of length L with the AP locations, where the real
%                  part is the horizontal position and the imaginary part
%                  is the vertical position
%UEpositions     = Vector of length K with UE positions, measured in the
%                  same way as APpositions
%distances       = Matrix with same dimension as gainOverNoisedB containing
%                  the distances in meter between APs and UEs
%
%This Matlab function was developed to generate simulation results to:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.
L = params.numGNB;
K_mmW = params.numUE;
K = params.numUE+params.numUE_sub6;
Lmax = params.Lmax;
N = params.num_antennas_per_gNB;
N_UE_FWA = params.N_UE_FWA;
N_UE_cell = params.N_UE_cell;
coverageRange = params.coverageRange;
coverageRange_sub6 = params.coverageRange_sub6;
% tau_p = params.tau_p;
ASD_varphi = params.ASD_varphi;
ASD_theta = params.ASD_theta;
%% Define simulation setup

%Set the seed number if it is specified other than zero
% if (nargin>9)&&(seed>0)
if (nargin>1)&&(seed>0)
    rng(seed)
end

% %Size of the coverage area (as a square with wrap-around)
% squareLength = 1000; %meter

%Communication bandwidth (Hz)
B = params.Band;
% B = params.scs_sub6;

%Noise figure (in dB)
noiseFigure = 7;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters for the model in (5.42)
alpha = 36.7;
constantTerm = -30.5;

%Standard deviation of the shadow fading in (5.43)
sigma_sf = 4;

%Decorrelation distance of the shadow fading in (5.43)
decorr = 9;

%Height difference between an AP and a UE (in meters)
% distanceVertical = 10;
distanceVertical = params.ht - params.hr;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Prepare to save results
gainOverNoisedB = zeros(L,K,nbrOfSetups);
R_gNB = zeros(N,N,L,K,nbrOfSetups);
R_ue_mmW = zeros(N_UE_FWA,N_UE_FWA,L,K_mmW,nbrOfSetups);
R_ue_sub6 = zeros(N_UE_cell,N_UE_cell,L,K-K_mmW,nbrOfSetups);
distances = zeros(L,K,nbrOfSetups);
pilotIndex = zeros(K,nbrOfSetups);
% D = zeros(L,K,nbrOfSetups);
D = ones(L,K,nbrOfSetups);
D_small = zeros(L,K,nbrOfSetups);

masterAPs = zeros(K,1); %the indices of master AP of each UE k 


%% Go through all setups
for n = 1:nbrOfSetups
    
    %Random AP locations with uniform distribution
    locationsBS = params.locationsBS;
    APpositions = locationsBS(:,1) + 1i*locationsBS(:,2);
    %Prepare to compute UE locations  
    UE_locations = params.UE_locations;
    UE_locations_sub6 = params.UE_locations_sub6;
    UEpositions = [UE_locations; UE_locations_sub6];
    UEpositions = UEpositions(:,1) + 1i*UEpositions(:,2);
    %Compute alternative AP locations by using wrap around
    % wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapHorizontal = repmat([-coverageRange_sub6 0 coverageRange_sub6],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

    %Prepare to store shadowing correlation matrix
    shadowCorrMatrix = sigma_sf^2*ones(K,K);
    shadowAPrealizations = zeros(K,L);    
    
    %Add UEs
    for k = 1:K
        
        %Generate a random UE location in the area
        % UEposition = (rand(1,1) + 1i*rand(1,1)) * squareLength;
        UEposition = UEpositions(k);        
        %Compute distances assuming that the APs are 10 m above the UEs
        [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - repmat(UEposition,size(APpositionsWrapped))),[],2);
        distances(:,k,n) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
        
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
            newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
            term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
            meanvalues = term1*shadowAPrealizations(1:k-1,:);
            stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
            
        else %If this is the first UE
            
            %Add the UE and begin to store shadow fading correlation values
            meanvalues = 0;
            stdvalue = sigma_sf;
            newcolumn = [];
            
        end
        
        %Generate the shadow fading realizations
        shadowing = meanvalues + stdvalue*randn(1,L);
        
        %Compute the channel gain divided by noise power
        gainOverNoisedB(:,k,n) = constantTerm - alpha*log10(distances(:,k,n)) + shadowing' - noiseVariancedBm;
        
        
        
        %Update shadowing correlation matrix and store realizations
        shadowCorrMatrix(1:k-1,k) = newcolumn;
        shadowCorrMatrix(k,1:k-1) = newcolumn';
        shadowAPrealizations(k,:) = shadowing;
        
        %Store the UE position
        UEpositions(k) = UEposition;
        
        
        %Determine the master AP for UE k by looking for AP with best
        %channel condition
%         [~,master] = max(gainOverNoisedB(:,k,n));
%         D(master,k,n) = 1;
%         masterAPs(k) = master;
%         
%         %Assign orthogonal pilots to the first tau_p UEs according to
%         %Algorithm 4.1
%         if k <= tau_p
%             
%             pilotIndex(k,n) = k;
%             
%         else %Assign pilot for remaining UEs
%             
%             %Compute received power to the master AP from each pilot
%             %according to Algorithm 4.1
%             pilotinterference = zeros(tau_p,1);
%             
%             for t = 1:tau_p
%                 
%                 pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1,n)==t,n)));
%                 
%             end
%             
%             %Find the pilot with the least receiver power according to
%             %Algorithm 4.1
%             [~,bestpilot] = min(pilotinterference);
%             pilotIndex(k,n) = bestpilot;
%             
%         end
        
        
        
        %Go through all APs
        for l = 1:L
            
            %Compute nominal angle between UE k and AP l
            angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
            angletoUE_theta = asin(distanceVertical/distances(l,k,n));  %elevation angle
            %Generate spatial correlation matrix using the local
            %scattering model in (2.18) and Gaussian angular distribution
            %by scaling the normalized matrices with the channel gain
            if nargin>12
                R_gNB(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*functionRlocalscattering_mod(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                if (k<=K_mmW)
                    R_ue_mmW(:,:,l,k,n) = functionRlocalscattering_mod(N_UE_FWA,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                else
                    R_ue_sub6(:,:,l,k-K_mmW,n) = functionRlocalscattering_mod(N_UE_cell,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
                end
            else
                R_gNB(:,:,l,k,n) = db2pow(gainOverNoisedB(l,k,n))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
                if (k<=K_mmW)
                    R_ue_mmW(:,:,l,k,n) = eye(N_UE_FWA);
                else
                    R_ue_sub6(:,:,l,k-K_mmW,n) = eye(N_UE_cell);
                end
            end
        end
        
    end
    
    
    %Each AP serves the UE with the strongest channel condition on each of
    %the pilots in the cell-free setup
%     for l = 1:L
%         
%         for t = 1:tau_p
%             
%             pilotUEs = find(t==pilotIndex(:,n));
%             [~,UEindex] = max(gainOverNoisedB(l,pilotUEs,n));
%             D(l,pilotUEs(UEindex),n) = 1;
%            
%         end
%         
%     end
    
    gainOverNoise = db2pow(gainOverNoisedB);
%     for l = 1:L
%         [gains, idxs] = sort(gainOverNoise(l,:), 'descend');
% %         for k = 1:K
% %             if ((sum(gains(1:k))/sum(gains))*100 > 90)
% %                 idxs_not_chosen = idxs((k+1):end);
% %                 break;
% %             end
% %         end
%         idxs_not_chosen = idxs((Kmax+1):end);
%         D(l,idxs_not_chosen,n) = 0;
%     end
%     for k = 1:K
%         if (sum(D(:,k,n)) == 0)
%             [~, idx] = max(gainOverNoise(:,k));
%             other_ue_idxs = setdiff(find(D(idx,:,n)),k);
%             [~,k_idx_idx] = min(gainOverNoise(idx,other_ue_idxs)); 
%             k_idx = other_ue_idxs(k_idx_idx);
%             D(idx,k,n) = 1;
%             D(idx,k_idx,n) = 0;
%         end
%     end
    for k = 1:K
        [gains, idxs] = sort(gainOverNoise(:,k), 'descend');
        if (k<=K_mmW)
            idxs_not_chosen = idxs(2:end);
        else
            idxs_not_chosen = idxs((Lmax+1):end);
        end
        D(idxs_not_chosen,k,n) = 0;
    end
    %Determine the AP serving each UE in the small-cell setup according to
    %(5.47) by considering only the APs from the set M_k for UE k, i.e.,
    %where D(:,k,n) is one.
    for k=1:K
        tempmat = -inf*ones(L,1);
        tempmat(D(:,k,n)==1,1) = gainOverNoisedB(D(:,k,n)==1,k,n);
        [~,servingAP] = max(tempmat);
        D_small(servingAP,k,n) = 1;        
    end    
end