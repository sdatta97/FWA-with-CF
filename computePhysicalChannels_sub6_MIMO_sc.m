function [phy_channel_sub6, phy_channel_sub6_est, phy_channel_mmW, phy_channel_mmW_est] = computePhysicalChannels_sub6_MIMO_sc(params)
K = params.numUE+params.numCPE;
K_FWA = params.numCPE;
Ntx = params.num_antennas_per_sc;
N_FWA = params.N_UE_FWA;
N_cell = params.N_UE_cell;
S = params.numSC;
BETA = params.BETA_sc;
R_gNB = params.R_sc;
R_cpe = params.R_cpe_sc;
R_ue = params.R_ue_sc;
phy_channel_mmW = zeros(S,K_FWA,Ntx,N_FWA);
phy_channel_sub6 = zeros(S,K-K_FWA,Ntx,N_cell);
for s = 1:S
    for k = 1:K_FWA
%         phy_channel_mmW (m,k,:,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,N_mmW) + 1i*randn(Ntx,N_mmW));        
        phy_channel_mmW (s,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,s,k))*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA))*sqrtm(R_cpe(:,:,s,k));        
    end
    for k = 1:K-K_FWA
%         phy_channel_sub6 (m,k,:,:) = sqrt(0.5*BETA(m,k+K_mmW))*(randn(Ntx,N_sub6) + 1i*randn(Ntx,N_sub6));        
        phy_channel_sub6 (s,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,s,k+K_FWA))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,s,k));        
    end 
end
if params.MOBILE
    for s = 1:S
        for k = 1:K-K_FWA
            phy_channel_sub6 (s,k,:,:) = params.mob_rho*phy_channel_sub6 (s,k,:,:) + reshape(sqrt(1 - params.mob_rho^2)*sqrt(0.5)*sqrtm(R_gNB(:,:,s,k+K_FWA,1))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,s,k,1)),[1,1,Ntx,N_cell]);        
        end 
    end
end
phy_channel_mmW_est = phy_channel_mmW;
phy_channel_sub6_est = phy_channel_sub6;
% phy_channel_sub6_est = zeros(M,Kd,Ntx);
% PHI1    = orth(rand(tau));   % generate an orthonormal matrix of dimension tau_p
% PHI     = zeros(size(PHI1));
% perm_vec  = repmat(randperm(tau),1,2);
% phi_index = perm_vec(1:K);
% for k = 1:K
%     PHI(:,k) = PHI1(:,phi_index(k));
% end
% for m = 1:M
%     W_tx = sqrt(0.5)*(randn(Ntx, tau)+1i*randn(Ntx,tau));
%     % W_rx = sqrt(0.5)*(randn(Nrx, tau)+1i*randn(Nrx,tau));
%     for k = 1:Kd
%         phy_channel_sub6 (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Ntx,1) + 1i*randn(Ntx,1));
%         % channel_bar_dl (m,k,:) = sqrt(0.5*BETA(m,k))*(randn(Nrx,1) + 1i*randn(Nrx,1));
%         c_dl (m,k) = sqrt(rho*tau)*BETA(m,k)/(1+rho*tau*BETA(m,k));
%         for n = 1:Ntx
%             phy_channel_sub6_est (m,k,n) = c_dl(m,k)*(sqrt(rho*tau)*phy_channel_sub6 (m,k,n)+ W_tx(n,:)*PHI(:,k));
%         end        
%     end
% end
end