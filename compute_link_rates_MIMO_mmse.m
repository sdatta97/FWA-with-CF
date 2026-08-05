function rate_dl = compute_link_rates_MIMO_mmse(params, channel_dl_FWA, channel_est_dl_FWA)
%FWA MU-MIMO rates on the shared FWA band: per-sector RZF/LP-MMSE
%precoding over the co-scheduled CPEs, with the fixed gamma_I residual
%interference suppression and the EVM hardware-impairment model. Called
%by the FWA phase of SimulationMain.m only (the cellular phase uses the
%compute_link_rates_OFDM* functions).
M_sectors = params.numGNB;
K_FWA = params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl_FWA,3);
N_CPE_FWA = size(channel_dl_FWA,4);
p_d = params.rho_tot; %total sector power. Constant-PSD convention: channels
%are normalized by FULL-band noise, so p_d*|h|^2 is the PSD ratio - exact
%on any subband; the bandwidth share enters only through the rate prelog
D_FWA = params.D_FWA;
SI_cancel_factor = params.SI_cancel_factor;
if params.HW_IMPAIRMENTS
    Kt     = params.Kt;
    Kr_FWA = params.Kr;
else
    Kt = 1; Kr_FWA = 1;
end

%Serving sector(s) of each CPE
Serv = cell(K_FWA,1);
for k = 1:K_FWA
    Serv{k} = find(D_FWA(:,k)==1);
end

%% effective-gain container and precoders
D_FWA_FWA = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
dl_mmse_precoder_FWA = zeros(size(channel_est_dl_FWA));
scaling_LP_mmse = zeros(M_sectors,K_FWA);
%MMSE precoding, structured like the OFDM beam matrices: the MMSE filter
%basis inv_matrix is identical for every user a sector precodes, so it
%is built once per sector and applied to each served user's channel.
%Beamforming operates on the ESTIMATED channels (imperfect CSI when
%params.IMPERFECT_CSI is set); physical reception below uses the true ones
for m = 1:M_sectors
    inv_matrix = eye(N_BS);
    for q = 1:K_FWA
        if ismember(m,Serv{q})
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])';
        end
    end
    for k = 1:K_FWA
        dl_mmse_precoder_FWA(m,k,:,:) = p_d*inv_matrix\(reshape(channel_est_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]));
        if ismember(m,Serv{k})
            scaling_LP_mmse(m,k) = norm(dl_mmse_precoder_FWA(m,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M_sectors
    for k = 1:K_FWA
        if ismember(m,Serv{k})
            dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])./sqrt(scaling_LP_mmse(m,k));
        end
    end
end
%Per-sector radiated-power constraint: eta_eq scales the precoders so
%sector m's TOTAL power over all served CPEs and streams is p_d
%(sum_q p_d*eta_eq(m,q)*trace(P_q P_q') = p_d by construction), the
%MU-MIMO counterpart of the cellular p_d/N_UE per-stream split
eta_eq = zeros(M_sectors,K_FWA);
for m = 1:M_sectors
    term = 0;
    for k = 1:K_FWA
        if ismember(m,Serv{k})
            term = term + trace(reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])');
        end
    end
    if (term > 0)
        eta_eq(m,:) = (1/term)*D_FWA(m,1:K_FWA);
    end
end
%Effective gains: the channels channel_dl_FWA already carry the sector
%antenna pattern of BS m towards each CPE (via the pattern-inclusive
%large-scale gains inside R_gNB), so both the desired links and every
%inter-sector interference term below are sectored
for k = 1:K_FWA
    for q = 1:K_FWA
        for m = 1:M_sectors
            if ismember(m,Serv{q})
                D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
            end
        end
    end
end

%% MIMO log-det rates
%Shannon capacity of user k's N_CPE x N_CPE effective MIMO channel
%G_k = D_FWA_FWA(k,k,:,:) under colored interference-plus-noise:
%  rate_k = BW*prelog*log2 det(I + Kt*Kr*G_k' R_k^{-1} G_k),
%R_k = sum_q!=k fac_q*Kt*Kr*G_kq G_kq' (spatial covariance of user q's
%precoded streams at k's antennas, fac = site-aware gamma_I) + EVM
%distortion (diagonal: per-antenna received power times (1-Kt*Kr),
%uncorrelated-distortion convention) + noise.
noise_FWA = abs(sqrt(0.5)*(randn(K_FWA,N_CPE_FWA) + 1j*randn(K_FWA,N_CPE_FWA))).^2;
%site of each CPE's serving sector (sectors are numbered site-major:
%params.locationsBS = kron(siteLocations, ones(S,1)))
siteAware = ~isfield(params,'gI_site_aware') || params.gI_site_aware; %0 = legacy blanket gamma_I (ablation)
S_site = params.sectors_per_site;
servSite = zeros(K_FWA,1);
for k = 1:K_FWA
    servSite(k) = ceil(Serv{k}(1)/S_site);
end
rate_dl = zeros(K_FWA,1);
for k = 1:K_FWA
    G = reshape(D_FWA_FWA(k,k,:,:),[N_CPE_FWA,N_CPE_FWA]); %rows = k's antennas, cols = k's streams
    R = diag(noise_FWA(k,:));                %interference-plus-noise covariance
    HIdiag = (1-Kr_FWA*Kt)*sum(abs(G).^2,2); %EVM distortion of the own signal
    for q = 1:K_FWA
        if (q~=k)
            %gamma_I models the CPE's DIRECTIONAL receive null and so
            %applies only to OTHER-SITE interferers: same-sector
            %streams are already suppressed by the L-MMSE precoder
            %(D_FWA_FWA carries the residual leakage), and co-site
            %other-sector interference arrives from the serving mast's
            %direction, inside the CPE's main receive lobe - neither
            %can be nulled again at the CPE
            if siteAware && servSite(q) == servSite(k), fac = 1; else, fac = SI_cancel_factor; end
            Gq = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]);
            R = R + Kr_FWA*Kt*fac*(Gq*Gq');
            HIdiag = HIdiag + (1-Kr_FWA*Kt)*fac*sum(abs(Gq).^2,2);
        end
    end
    R = R + diag(HIdiag);
    rate_dl(k) = BW*TAU_FAC*logdet2(Kr_FWA*Kt*(G'*(R\G)));
end
end
