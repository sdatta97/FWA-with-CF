function [phy_channel_ue, phy_channel_ue_est, phy_channel_FWA, phy_channel_FWA_est, phy_channel_interFWA, phy_channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params)
M_sectors = params.numGNB; %number of sector-cell BS entries (M_sectors in the main scripts)
K = M_sectors*params.numUE+params.numCPE;
K_FWA = params.numCPE;
Ntx = params.num_antennas_per_gNB;
N_FWA = params.N_UE_FWA;
N_cell = params.N_UE_cell;
%NOTE: R_gNB / R_cpe / R_interue are CHOLESKY FACTORS L (R = L*L') of
%the spatial correlation matrices, so channels are drawn as L*Z*L2' with
%no per-realization matrix square roots; R(1,1) = |L(1,1)|^2 recovers the
%per-element channel variance beta for the LMMSE-CSI block
R_gNB = params.R_gNB;
R_cpe = params.R_cpe;
R_interue = params.R_interue;
R_ue = params.R_ue;
phy_channel_FWA = zeros(M_sectors,K_FWA,Ntx,N_FWA);
phy_channel_ue = zeros(M_sectors,K-K_FWA,Ntx,N_cell);
phy_channel_interFWA = zeros(K_FWA,K,N_FWA,N_FWA);
for m = 1:M_sectors
    for k = 1:K_FWA
        phy_channel_FWA (m,k,:,:) = sqrt(0.5)*R_gNB(:,:,m,k)*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA))*R_cpe(:,:,m,k)';
    end
    for k = 1:K-K_FWA
        phy_channel_ue (m,k,:,:) = sqrt(0.5)*R_gNB(:,:,m,k+K_FWA)*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*R_ue(:,:,m,k)';
    end
end
if params.MOBILE
    %Channel aging factor per cellular UE (indoor and outdoor UEs move at
    %different velocities); a scalar is replicated for backward compatibility
    mob_rho = params.mob_rho;
    if isscalar(mob_rho)
        mob_rho = mob_rho*ones(K-K_FWA,1);
    end
    for m = 1:M_sectors
        for k = 1:K-K_FWA
            phy_channel_ue (m,k,:,:) = mob_rho(k)*phy_channel_ue (m,k,:,:) + reshape(sqrt(1 - mob_rho(k)^2)*sqrt(0.5)*R_gNB(:,:,m,k+K_FWA)*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell))*R_ue(:,:,m,k)', [1,1,Ntx,N_cell]);
        end
    end
end
for q = 1:K_FWA
    for k = 1:K
        if (q~=k)
            phy_channel_interFWA(q,k,:,:) = sqrt(0.5)*R_interue(:,:,q,k)*(randn(N_FWA,N_FWA) + 1i*randn(N_FWA,N_FWA));
        end
    end
end
if params.IMPERFECT_CSI
    %Imperfect CSI: each link is estimated from CSI-RS observations
    %(TS 38.214) received in white noise, followed by per-link LMMSE
    %estimation. In the noise-normalized units used throughout, the
    %pilot observation of a channel entry h ~ CN(0,beta) is
    %y = h + n/sqrt(p_csi) with n ~ CN(0,1), so the LMMSE estimate is
    %h_est = (gamma/(1+gamma))*y with pilot SNR gamma = p_csi*beta.
    %CSI-RS EPRE follows the PDSCH EPRE plus the standardised
    %powerControlOffset (TS 38.214, range -8..15 dB, default 0 dB).
    p_csi = params.rho_tot*10^(0.1*params.csi_rs_offset_dB);
    phy_channel_FWA_est = zeros(size(phy_channel_FWA));
    phy_channel_ue_est = zeros(size(phy_channel_ue));
    phy_channel_interFWA_est = zeros(size(phy_channel_interFWA));
    for m = 1:M_sectors
        for k = 1:K_FWA
            beta = abs(R_gNB(1,1,m,k))^2; %per-element variance from the factor diagonal
            c = p_csi*beta/(p_csi*beta + 1);
            noise = sqrt(0.5/p_csi)*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA));
            phy_channel_FWA_est(m,k,:,:) = c*(reshape(phy_channel_FWA(m,k,:,:),[Ntx,N_FWA]) + noise);
        end
        for k = 1:K-K_FWA
            beta = abs(R_gNB(1,1,m,k+K_FWA))^2; %factor diagonal
            c = p_csi*beta/(p_csi*beta + 1);
            noise = sqrt(0.5/p_csi)*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell));
            phy_channel_ue_est(m,k,:,:) = c*(reshape(phy_channel_ue(m,k,:,:),[Ntx,N_cell]) + noise);
        end
    end
    %repeater-to-user links: sounded through the network-controlled
    %repeater with the same normalized pilot power
    for q = 1:K_FWA
        for k = 1:K
            if (q~=k)
                beta = abs(R_interue(1,1,q,k))^2; %factor diagonal
                c = p_csi*beta/(p_csi*beta + 1);
                noise = sqrt(0.5/p_csi)*(randn(N_FWA,N_FWA) + 1i*randn(N_FWA,N_FWA));
                phy_channel_interFWA_est(q,k,:,:) = c*(reshape(phy_channel_interFWA(q,k,:,:),[N_FWA,N_FWA]) + noise);
            end
        end
    end
else
    %perfect CSI: the estimates equal the true channels
    phy_channel_FWA_est = phy_channel_FWA;
    phy_channel_ue_est = phy_channel_ue;
    phy_channel_interFWA_est = phy_channel_interFWA;
end
end
