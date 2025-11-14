function [phy_channel_sub6, phy_channel_sub6_est, phy_channel_mmW, phy_channel_mmW_est] = computePhysicalChannels_sub6_MIMO(params)
K = params.numUE+params.numCPE;
K_FWA = params.numCPE;
Ntx = params.num_antennas_per_gNB;
N_FWA = params.N_UE_FWA;
N_cell = params.N_UE_cell;
M = params.numGNB;
BETA = params.BETA;
R_gNB = params.R_gNB;
R_cpe = params.R_cpe;
R_ue = params.R_ue;
phy_channel_mmW = zeros(M,K_FWA,Ntx,N_FWA);
phy_channel_sub6 = zeros(M,K-K_FWA,Ntx,N_cell);
for m = 1:M
    for k = 1:K_FWA
        phy_channel_mmW (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k,1))*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA))*sqrtm(R_cpe(:,:,m,k,1));        
    end
    for k = 1:K-K_FWA
        phy_channel_sub6 (m,k,:,:) = sqrt(0.5)*sqrtm(R_gNB(:,:,m,k+K_FWA,1))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,m,k,1));        
    end 
end
if params.MOBILE
    for m = 1:M
        for k = 1:K-K_FWA
            phy_channel_sub6 (m,k,:,:) = params.mob_rho*phy_channel_sub6 (m,k,:,:) + reshape(sqrt(1 - params.mob_rho^2)*sqrt(0.5)*sqrtm(R_gNB(:,:,m,k+K_FWA,1))*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*sqrtm(R_ue(:,:,m,k,1)), [1,1,Ntx,N_cell]);        
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