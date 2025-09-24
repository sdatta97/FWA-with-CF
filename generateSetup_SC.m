function [gainOverNoisedB,R_gNB,R_ue_mmW,R_ue_sub6,D,D_small,APpositions,distances] = generateSetup_SC(params,seed,pilotIndex,D,D_small,UEpositions)
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
%                  element (l,k) is the channel gain (normalized by the
%                  noise variance) between AP l and UE k in setup n
%R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                  where (:,:,l,k) is the spatial correlation matrix
%                  between AP l and UE k in setup n, normalized by noise
%pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                  pilots assigned to the UEs
%D               = DCC matrix with dimension L x K x nbrOfSetups where (l,k)
%                  is one if AP l serves UE k in setup n and zero otherwise
%                  for cell-free setup
%D_small         = DCC matrix with dimension L x K x nbrOfSetups where (l,k)
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
L = size(params.locationsSC,1);
K_FWA = params.numCPE;
K = params.numUE+params.numCPE;
Lmax = params.Lmax;
N = params.num_antennas_per_sc;
N_UE_FWA = params.N_UE_FWA;
N_UE_cell = params.N_UE_cell;
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
distanceVertical = params.ht_sc - params.hr;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance


%Prepare to save results
gainOverNoisedB = zeros(L,K);
R_gNB = zeros(N,N,L,K);
R_ue_mmW = zeros(N_UE_FWA,N_UE_FWA,L,K_FWA);
R_ue_sub6 = zeros(N_UE_cell,N_UE_cell,L,K-K_FWA);
distances = zeros(L,K);
D = ones(L,K);
D_small = zeros(L,K);
masterAPs = zeros(K,1); %the indices of master AP of each UE k 


%% Go through all setups
    
%Random AP locations with uniform distribution
locationsSC = params.locationsSC;
APpositions = locationsSC(:,1) + 1i*locationsSC(:,2);
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
    distances(:,k) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
    
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
    gainOverNoisedB(:,k) = constantTerm - alpha*log10(distances(:,k)) + shadowing' - noiseVariancedBm;
    
    
    
    %Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:k-1,k) = newcolumn;
    shadowCorrMatrix(k,1:k-1) = newcolumn';
    shadowAPrealizations(k,:) = shadowing;
    
    
    %Determine the master AP for UE k by looking for AP with best
    %channel condition
%         [~,master] = max(gainOverNoisedB(:,k));
%         D(master,k) = 1;
%         masterAPs(k) = master;
%         
%         %Assign orthogonal pilots to the first tau_p UEs according to
%         %Algorithm 4.1
%         if k <= tau_p
%             
%             pilotIndex(k) = k;
%             
%         else %Assign pilot for remaining UEs
%             
%             %Compute received power to the master AP from each pilot
%             %according to Algorithm 4.1
%             pilotinterference = zeros(tau_p,1);
%             
%             for t = 1:tau_p
%                 
%                 pilotinterference(t) = sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1)==t)));
%                 
%             end
%             
%             %Find the pilot with the least receiver power according to
%             %Algorithm 4.1
%             [~,bestpilot] = min(pilotinterference);
%             pilotIndex(k) = bestpilot;
%             
%         end
    
    
    
    %Go through all APs
    for l = 1:L
        
        %Compute nominal angle between UE k and AP l
        angletoUE_varphi = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); %azimuth angle
        angletoUE_theta = asin(distanceVertical/distances(l,k));  %elevation angle
        %Generate spatial correlation matrix using the local
        %scattering model in (2.18) and Gaussian angular distribution
        %by scaling the normalized matrices with the channel gain
        if nargin>12
            R_gNB(:,:,l,k) = db2pow(gainOverNoisedB(l,k))*functionRlocalscattering_mod(N,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            if (k<=K_FWA)
                R_ue_mmW(:,:,l,k) = functionRlocalscattering_mod(N_UE_FWA,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            else
                R_ue_sub6(:,:,l,k-K_FWA) = functionRlocalscattering_mod(N_UE_cell,angletoUE_varphi,angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
            end
        else
            R_gNB(:,:,l,k) = db2pow(gainOverNoisedB(l,k))*eye(N);  %If angular standard deviations are not specified, set i.i.d. fading
            if (k<=K_FWA)
                R_ue_mmW(:,:,l,k) = eye(N_UE_FWA);
            else
                R_ue_sub6(:,:,l,k-K_FWA) = eye(N_UE_cell);
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
%             pilotUEs = find(t==pilotIndex(:));
%             [~,UEindex] = max(gainOverNoisedB(l,pilotUEs));
%             D(l,pilotUEs(UEindex)) = 1;
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
%         D(l,idxs_not_chosen) = 0;
%     end
%     for k = 1:K
%         if (sum(D(:,k)) == 0)
%             [~, idx] = max(gainOverNoise(:,k));
%             other_ue_idxs = setdiff(find(D(idx,:)),k);
%             [~,k_idx_idx] = min(gainOverNoise(idx,other_ue_idxs)); 
%             k_idx = other_ue_idxs(k_idx_idx);
%             D(idx,k) = 1;
%             D(idx,k_idx) = 0;
%         end
%     end
for k = 1:K
    [gains, idxs] = sort(gainOverNoise(:,k), 'descend');
    idxs_not_chosen = idxs((Lmax+1):end);
    D(idxs_not_chosen,k) = 0;
end
%Determine the AP serving each UE in the small-cell setup according to
%(5.47) by considering only the APs from the set M_k for UE k, i.e.,
%where D(:,k) is one.
for k=1:K
    tempmat = -inf*ones(L,1);
    tempmat(D(:,k)==1,1) = gainOverNoisedB(D(:,k)==1,k);
    [~,servingAP] = max(tempmat);
    D_small(servingAP,k) = 1;        
end    
end