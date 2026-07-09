function [phy_channel_ue, phy_channel_ue_est, phy_channel_FWA, phy_channel_FWA_est, phy_channel_interFWA, phy_channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params)
M = params.numGNB; %number of sector-cell BS entries (M_sectors in the main scripts)
K = M*params.numUE+params.numCPE;
K_FWA = params.numCPE;
Ntx = params.num_antennas_per_gNB;
N_FWA = params.N_UE_FWA;
N_cell = params.N_UE_cell;
R_gNB = params.R_gNB;
R_cpe = params.R_cpe;
R_interue = params.R_interue;
R_ue = params.R_ue;
phy_channel_FWA = zeros(M,K_FWA,Ntx,N_FWA);
phy_channel_ue = zeros(M,K-K_FWA,Ntx,N_cell);
phy_channel_interFWA = zeros(K_FWA,K,N_FWA,N_FWA);
for m = 1:M
    for k = 1:K_FWA
        phy_channel_FWA (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k))*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA))*sqrtm(R_cpe(:,:,m,k));     
    end
    for k = 1:K-K_FWA
        phy_channel_ue (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k+K_FWA))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,m,k));
    end 
end
if params.MOBILE
    %Channel aging factor per cellular UE (indoor and outdoor UEs move at
    %different velocities); a scalar is replicated for backward compatibility
    mob_rho = params.mob_rho;
    if isscalar(mob_rho)
        mob_rho = mob_rho*ones(K-K_FWA,1);
    end
    for m = 1:M
        for k = 1:K-K_FWA
            phy_channel_ue (m,k,:,:) = mob_rho(k)*phy_channel_ue (m,k,:,:) + reshape(sqrt(1 - mob_rho(k)^2)*sqrt(0.5)*sqrtm(R_gNB(:,:,m,k+K_FWA))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,m,k)), [1,1,Ntx,N_cell]);
        end
    end
end
for q = 1:K_FWA
    for k = 1:K
        if (q~=k)
            phy_channel_interFWA(q,k,:,:) = sqrt(0.5)*sqrtm(R_interue(:,:,q,k))*(randn(N_FWA,N_FWA) + 1i*randn(N_FWA,N_FWA));
        end
    end
end
phy_channel_FWA_est = phy_channel_FWA;
phy_channel_ue_est = phy_channel_ue;
phy_channel_interFWA_est = phy_channel_interFWA;
% perm_vec  = repmat(randperm(tau),1,2);
% phi_index = perm_vec(1:K);
%     W_tx = sqrt(0.5)*(randn(Ntx, tau)+1i*randn(Ntx,tau));
%     % W_rx = sqrt(0.5)*(randn(Nrx, tau)+1i*randn(Nrx,tau));
%         % channel_bar_dl (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Nrx,1) + 1i*randn(Nrx,1));
%         c_dl (m,k) = sqrt(rho*tau)*BETA(m,k)/(1+rho*tau*BETA(m,k));
end