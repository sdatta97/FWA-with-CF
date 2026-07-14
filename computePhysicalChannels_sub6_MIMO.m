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
%TRUE CSI aging (Truong & Heath, https://arxiv.org/abs/1305.6151): the
%drawn phy_channel_ue is the SOUNDED channel the gNB estimates from; the
%channel actually seen at data transmission has evolved by the per-UE
%Jakes factor rho = J0(2*pi*f_D*T_loop/2) set in SimulationMain. The
%sounded copy is kept for the estimation block below, so the beams are
%computed from STALE observations while reception uses the evolved truth.
phy_channel_ue_snd = phy_channel_ue;
if params.MOBILE
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
    contam = isfield(params,'fwa_pilot_contam') && params.fwa_pilot_contam ...
        && isfield(params,'pilot_fwa') && numel(params.pilot_fwa) == K_FWA;
    if contam
        %PILOT-CONTAMINATED SRS estimation of the BS-to-CPE links (TDD
        %reciprocity, TS 38.214): every gNB sector correlates its uplink
        %receive signal with each (sequence root, cyclic shift); CPEs of
        %the same site sharing a pilot land in ONE observation Y, so
        %their estimates are collinear scaled copies (the contamination
        %that poisons MU-MIMO precoding); the other site's same-comb-
        %offset CPEs leak into Y at the cross-root correlation power.
        %LMMSE per element: c_k = p*beta_k / (p*(sum of co-pilot betas +
        %leak*sum of cross-site co-offset betas) + 1). Pilot plan and
        %standards references in SimulationMain.
        pil = params.pilot_fwa(:);
        site_cpe = params.site_of_cpe(:);
        leak = 10^(0.1*params.srs_root_leak_dB);
        n_shift = 12; %cyclic shifts per comb offset (comb-4, TS 38.211)
        comb_of = ceil(pil/n_shift);
        beta_fwa = zeros(M_sectors,K_FWA);
        for m = 1:M_sectors
            for k = 1:K_FWA
                beta_fwa(m,k) = abs(R_gNB(1,1,m,k))^2;
            end
        end
        M_sites_p = max(site_cpe);
        for m = 1:M_sectors
            for s = 1:M_sites_p
                for p = unique(pil(site_cpe == s))'
                    members = find(site_cpe == s & pil == p);
                    xleak = find(site_cpe ~= s & comb_of == ceil(p/n_shift));
                    Y = sqrt(0.5/p_csi)*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA));
                    for q = members'
                        Y = Y + reshape(phy_channel_FWA(m,q,:,:),[Ntx,N_FWA]);
                    end
                    for q = xleak'
                        Y = Y + sqrt(leak)*reshape(phy_channel_FWA(m,q,:,:),[Ntx,N_FWA]);
                    end
                    den = p_csi*(sum(beta_fwa(m,members)) + leak*sum(beta_fwa(m,xleak))) + 1;
                    for k = members'
                        phy_channel_FWA_est(m,k,:,:) = (p_csi*beta_fwa(m,k)/den)*Y;
                    end
                end
            end
        end
    end
    for m = 1:M_sectors
        if ~contam
        for k = 1:K_FWA
            beta = abs(R_gNB(1,1,m,k))^2; %per-element variance from the factor diagonal
            c = p_csi*beta/(p_csi*beta + 1);
            noise = sqrt(0.5/p_csi)*(randn(Ntx,N_FWA) + 1i*randn(Ntx,N_FWA));
            phy_channel_FWA_est(m,k,:,:) = c*(reshape(phy_channel_FWA(m,k,:,:),[Ntx,N_FWA]) + noise);
        end
        end
        for k = 1:K-K_FWA
            beta = abs(R_gNB(1,1,m,k+K_FWA))^2; %factor diagonal
            c = p_csi*beta/(p_csi*beta + 1);
            noise = sqrt(0.5/p_csi)*(randn(Ntx,N_cell) + 1i*randn(Ntx,N_cell));
            %estimate from the SOUNDED channel (stale by the aging above)
            phy_channel_ue_est(m,k,:,:) = c*(reshape(phy_channel_ue_snd(m,k,:,:),[Ntx,N_cell]) + noise);
        end
    end
    %repeater-to-user links: sounded through the network-controlled
    %repeater with the same normalized pilot power. DESIGN ASSUMPTION
    %(beyond Rel-18): the SMART NCR measures its service channel locally
    %at the fading rate, so these estimates carry no staleness even
    %toward fast UEs - unlike the gNB links above, which age with the
    %per-class CSI loops. This freshness is what lets the repeater path
    %bypass the mobility penalty of the direct link (a deliberate,
    %standard-defying capability of the modeled repeater; Rel-18 side
    %control alone is semi-static, TR 38.867)
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
