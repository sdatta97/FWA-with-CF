function rate_dl = compute_link_rates_MIMO_mmse_wi_repeater(params,channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est)
%MU-MIMO MMSE downlink rates (FWA phase) WITH FWA-side NCR assistance,
%selected by the params.FWA_REPEAT flag in SimulationMain.m (the
%repeater-free variant lives in compute_link_rates_MIMO_mmse.m).
%When params.fwa_rep_pool is
%nonempty and the inter-CPE channels are supplied, the CPE-hosted NCRs
%ALSO assist the FWA band: the weakest params.rep_assist_frac fraction
%of each sector's CPEs may attach to the SAME enabled repeater pool as
%the cellular side, under the same donor-sector lock, ZF donor
%combining, output-power cap, and per-victim benefit gate
%(params.ncr_benefit_gate) as compute_link_rates_OFDM_wi_repeater.m.
%Differences forced by the MU-MIMO FWA phase (all CPEs are spatially
%multiplexed over the whole band, unlike per-UE OFDMA subbands):
% - an NCR with any attachment forwards the FULL FWA band, so the
%   33 dBm cap applies to the FWA-band donor input power
%   params.rep_P_in_FWA_mW (constant-PSD share, supplied by the caller);
% - a CPE never attaches to the NCR mounted on itself;
% - forwarded components on NON-assisted CPEs' channels are neglected
%   (mirror of the cellular OFF-variant approximation), and the AF-chain
%   distortion of the forwarded part enters as a sqrt(Kt_rep*Kr_rep)
%   coherent-amplitude loss (EVM-level; the residual distortion power is
%   neglected).
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
N_BS = size(channel_dl,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot;
if params.HW_IMPAIRMENTS
    Kt = params.Kt; Kr_FWA = params.Kr; Kr_cell = params.Kr;
    K_rep_hw = params.Kt_rep*params.Kr_rep;
else
    Kt = 1; Kr_FWA = 1; Kr_cell = 1; K_rep_hw = 1;
end
%noise draws shared by the NCR-OFF/ON evaluations; nothing above them
%consumes randomness, so the stream position matches the previous
%implementation exactly
noise_FWA = abs(sqrt(0.5)*(randn(K_FWA,N_CPE_FWA) + 1j*randn(K_FWA,N_CPE_FWA))).^2;
noise_Cell = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
%% FWA-side NCR assistance: composite serving channels for assisted CPEs
assisted = false(K_FWA,1);
if nargin >= 7 && isfield(params,'fwa_rep_pool') && ~isempty(params.fwa_rep_pool) && K_FWA > 0
    pool = params.fwa_rep_pool(:)';
    D_FWA_a = params.D_FWA;
    donor_of = zeros(K_FWA,1);
    for q = 1:K_FWA
        donor_of(q) = find(D_FWA_a(:,q)==1,1);
    end
    %eligibility mirrors the cellular side: the weakest rep_assist_frac
    %fraction of each sector's CPEs by serving-link large-scale gain
    beta_cpe = params.gainOverNoise_lin(:,1:K_FWA);
    elig = false(K_FWA,1);
    for m = 1:M_sectors
        served = find(donor_of==m)';
        served = served(~ismember(served, params.set_repeat)); %not scheduled
            %on this subband (hybrid split) -> not eligible for assistance
        if isempty(served), continue, end
        [~,ord] = sort(beta_cpe(m,served),'ascend');
        n_assist = ceil(params.rep_assist_frac*numel(served));
        elig(served(ord(1:n_assist))) = true;
    end
    %attachment: donor-sector-locked, at most num_repeater_per_cpe hosts,
    %and never the NCR mounted on the CPE itself
    BETA_interCPE = params.BETA_interUE(1:K_FWA,1:K_FWA);
    RepF = cell(K_FWA,1);
    n_att = zeros(K_FWA,1);
    for k = 1:K_FWA
        if elig(k)
            v = -Inf(K_FWA,1);
            for q = pool
                if q ~= k && donor_of(q) == donor_of(k)
                    v(q) = beta_cpe(donor_of(q),q)*params.rep_donor_gain(donor_of(q),q)*BETA_interCPE(q,k);
                end
            end
            [vmax,hosts] = maxk(v,params.num_repeater_per_cpe);
            RepF{k} = hosts(isfinite(vmax));
            n_att(RepF{k}) = n_att(RepF{k}) + 1;
        end
    end
    assisted = ~cellfun(@isempty,RepF);
end
if any(assisted)
    if isfield(params,'ncr_rank')
        ncr_rank = min([params.ncr_rank, N_CPE_FWA]);
    else
        ncr_rank = 1;
    end
    G_max_lin = 10^(0.1*params.repeat_gain);
    P_max_mW = 10^(0.1*params.rep_max_pow_dBm);
    rep_gain_amp = zeros(K_FWA,1);
    for q = find(n_att>0)'
        %full-FWA-band forwarding (forwarded fraction 1, see header)
        rep_gain_amp(q) = sqrt(min(G_max_lin, P_max_mW/params.rep_P_in_FWA_mW(q)));
    end
    rep_donor_amp = sqrt(params.rep_donor_gain);
    %per-repeater donor receive beams: ZF/max-SINR generalized-eigen
    %combining, identical to the cellular side
    W_r_rep = cell(K_FWA,1);
    for q = find(n_att>0)'
        m_donor = donor_of(q);
        Hd_est = reshape(channel_est_dl_FWA(m_donor,q,:,:),[N_BS,N_CPE_FWA]);
        if ~isfield(params,'ncr_zf_donor') || params.ncr_zf_donor
            Bq = eye(N_CPE_FWA);
            for mm = setdiff(1:M_sectors,m_donor)
                Hi = reshape(channel_est_dl_FWA(mm,q,:,:),[N_BS,N_CPE_FWA]);
                Bq = Bq + (p_d/N_BS)*params.rep_donor_gain(mm,q)*(Hi'*Hi);
            end
            [Vg,Dg] = eig(Hd_est'*Hd_est, Bq, 'vector');
            [~,ordg] = sort(real(Dg),'descend');
            Vr = Vg(:,ordg);
            Vr = Vr./vecnorm(Vr);
        else
            [~,~,Vr] = svd(Hd_est,'econ');
        end
        W_r_rep{q} = Vr(:,1:min(ncr_rank,size(Vr,2)));
    end
    ch_dl_FWA_on = channel_dl_FWA;
    ch_est_FWA_on = channel_est_dl_FWA;
    for k = find(assisted)'
        for kk = 1:numel(RepF{k})
            q = RepF{k}(kk);
            Vr = W_r_rep{q};
            Hs_est = reshape(channel_interFWA_est(q,k,:,1:N_CPE_FWA),[N_CPE_FWA,N_CPE_FWA]);
            [~,~,Vt] = svd(Hs_est.','econ');
            for chain = 1:ncr_rank
                w_r = Vr(:,chain); w_t = Vt(:,chain);
                amp_chain = rep_gain_amp(q)/sqrt(ncr_rank);
                for m = 1:M_sectors
                    g_vec = reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,q)*amp_chain;
                    g_vec_est = reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,q)*amp_chain;
                    for n = 1:N_CPE_FWA
                        s_n = reshape(channel_interFWA(q,k,:,n),[1,N_CPE_FWA])*w_t;
                        s_n_est = reshape(channel_interFWA_est(q,k,:,n),[1,N_CPE_FWA])*w_t;
                        ch_dl_FWA_on(m,k,:,n) = reshape(ch_dl_FWA_on(m,k,:,n),[N_BS,1]) + sqrt(K_rep_hw)*g_vec*s_n;
                        ch_est_FWA_on(m,k,:,n) = reshape(ch_est_FWA_on(m,k,:,n),[N_BS,1]) + g_vec_est*s_n_est;
                    end
                end
            end
        end
    end
end
%NCR-OFF rates (the previous implementation, verbatim, in rates_core_mmse)
rate_dl = rates_core_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, noise_FWA, noise_Cell, Kt, Kr_FWA, Kr_cell);
if any(assisted)
    %NCR-ON rates with the assisted CPEs' composite channels; the per-CPE
    %benefit gate keeps the better of ON/OFF (side control), exactly as
    %on the cellular side
    rate_on = rates_core_mmse(params, channel_dl, channel_est_dl, ch_dl_FWA_on, ch_est_FWA_on, noise_FWA, noise_Cell, Kt, Kr_FWA, Kr_cell);
    gate_on = ~isfield(params,'ncr_benefit_gate') || params.ncr_benefit_gate;
    for k = find(assisted)'
        if rate_on(k) > rate_dl(k) || ~gate_on
            rate_dl(k) = rate_on(k);
        end
    end
end
end

function rate_dl = rates_core_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, noise_FWA, noise_Cell, Kt, Kr_FWA, Kr_cell)
%the pre-assistance compute_link_rates_MIMO_mmse computation, verbatim,
%with the channels and shared noise draws passed in
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; %total sector power. Constant-PSD convention: channels
%are normalized by FULL-band noise, so p_d*|h|^2 is the PSD ratio - exact
%on any subband; the bandwidth share enters only through the rate prelog
D_FWA = params.D_FWA;
set_repeat = params.set_repeat;
SI_cancel_factor = params.SI_cancel_factor;
%Serving sector(s) of each user
Serv = cell(K,1);
for k = 1:K
    Serv{k} = find(D_FWA(:,k)==1);
end

%% initialization of c
D_FWA_FWA = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
D_FWA_Cell = zeros(K_FWA,K-K_FWA,N_CPE_FWA,N_UE);
D_Cell_FWA = zeros(K-K_FWA,K_FWA,N_UE,N_CPE_FWA);
D_Cell_Cell = zeros(K-K_FWA,K-K_FWA,N_UE,N_UE);
dl_mmse_precoder_FWA = zeros(size(channel_est_dl_FWA));
dl_mmse_precoder = zeros(size(channel_est_dl));
scaling_LP_mmse = zeros(M_sectors,K);
%MMSE precoding, structured like the OFDM beam matrices: the MMSE filter
%basis inv_matrix is identical for every user a sector precodes, so it
%is built once per sector and applied to each served user's channel.
%Beamforming operates on the ESTIMATED channels (imperfect CSI when
%params.IMPERFECT_CSI is set); physical reception below uses the true ones
for m = 1:M_sectors
    inv_matrix = eye(N_BS);
    for q = 1:K_FWA
        if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])';
        end
    end
    for q = 1:K-K_FWA
        if ismember(m,Serv{q+K_FWA})
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])';
        end
    end
    for k = 1:K_FWA
        if ~ismember(k,set_repeat)
            dl_mmse_precoder_FWA(m,k,:,:) = p_d*inv_matrix\(reshape(channel_est_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]));
            if ismember(m,Serv{k})
                scaling_LP_mmse(m,k) = norm(dl_mmse_precoder_FWA(m,k,:,:),'fro')^2;
            end
        end
    end
    for k = 1:K-K_FWA
        dl_mmse_precoder(m,k,:,:) = p_d*inv_matrix\(reshape(channel_est_dl(m,k,:,:),[N_BS,N_UE]));
        if ismember(m,Serv{k+K_FWA})
            scaling_LP_mmse(m,k+K_FWA) = norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M_sectors
    for k = 1:K_FWA
        if ismember(m,Serv{k}) && ~ismember(k,set_repeat)
            dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])./sqrt(scaling_LP_mmse(m,k));
        end
    end
    for k = 1:K-K_FWA
        if ismember(m,Serv{k+K_FWA})
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE])./sqrt(scaling_LP_mmse(m,k+K_FWA));
        end
    end
end
eta_eq = zeros(M_sectors,K);
for m = 1:M_sectors
    term = 0;
    for k = 1:K
        if ismember(m,Serv{k})
            if (k<=K_FWA)
                if ~ismember(k,set_repeat)
                    term = term + trace(reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])');
                end
            else
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])*reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])');
            end
        end
    end
    if (term > 0)
        setA = 1:K;
        eta_eq(m,setA(~ismember(setA,set_repeat))) = (1/term)*D_FWA(m,setA(~ismember(setA,set_repeat)));
    end
end
%Effective gains: the channels channel_dl_FWA/channel_dl already carry
%the sector antenna pattern of BS m towards each user (via the
%pattern-inclusive large-scale gains inside R_gNB), so both the desired
%links and every inter-sector interference term below are sectored
for k = 1:K_FWA
    if ~ismember(k,set_repeat)
        for q = 1:K_FWA
            for m = 1:M_sectors
                if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
                    D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                end
            end
        end
        for q = 1:K-K_FWA
            for m = 1:M_sectors
                if ismember(m,Serv{q+K_FWA})
                    D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
                end
            end
        end
    end
end
for k = 1:K-K_FWA
    for q = 1:K_FWA
        for m = 1:M_sectors
            if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
                D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
            end
        end
    end
    for q = 1:K-K_FWA
        for m = 1:M_sectors
            if ismember(m,Serv{q+K_FWA})
                D_Cell_Cell(k,q,:,:) = reshape(D_Cell_Cell(k,q,:,:),[N_UE,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
            end
        end
    end
end
DS_FWA = zeros(K_FWA,N_CPE_FWA);
MSI_FWA = zeros(K_FWA,N_CPE_FWA);
MUI_FWA = zeros(K_FWA,N_CPE_FWA);
HI_FWA = zeros(K_FWA,N_CPE_FWA);
DS_Cell = zeros(K-K_FWA,N_UE);
MSI_Cell = zeros(K-K_FWA,N_UE);
MUI_Cell = zeros(K-K_FWA,N_UE);
HI_Cell = zeros(K-K_FWA,N_UE);

snr_num_FWA = zeros(K_FWA,N_CPE_FWA);
snr_den_FWA = zeros(K_FWA,N_CPE_FWA);
snr_num_Cell = zeros(K-K_FWA,N_UE);
snr_den_Cell = zeros(K-K_FWA,N_UE);
rate_dl = zeros(K,1);
for k = 1:K_FWA
    if ~ismember(k,set_repeat)
        for n = 1:N_CPE_FWA
            DS_FWA(k,n)  = Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,n)))^2;
            HI_FWA(k,n)  = (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,n)))^2;
            for nn = 1:N_CPE_FWA
                if (abs(D_FWA_FWA(k,k,nn,nn))<abs(D_FWA_FWA(k,k,n,n)))
                    MSI_FWA(k,n) = MSI_FWA(k,n) + Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                    HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                end
            end
            for q = 1:K_FWA
                if (q~=k) && ~ismember(q,set_repeat)
                    MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*SI_cancel_factor*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                    HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*SI_cancel_factor*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                end
            end
            for q = 1:K-K_FWA
                MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
                HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
            end
            snr_num_FWA(k,n) = DS_FWA(k,n);
            snr_den_FWA(k,n) = MSI_FWA(k,n) + MUI_FWA(k,n) + HI_FWA(k,n) + noise_FWA(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_FWA(k,n)/snr_den_FWA(k,n));
        end
    end
end
for k = 1:K-K_FWA
    for n = 1:N_UE
        DS_Cell(k,n)  = Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,n)))^2;
        HI_Cell(k,n)  = (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,n)))^2;
        for nn = 1:N_UE
            if (abs(D_Cell_Cell(k,k,nn,nn))<abs(D_Cell_Cell(k,k,n,n)))
                MSI_Cell(k,n) = MSI_Cell(k,n) + Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,nn)))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_FWA
            if ~ismember(q,set_repeat)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
            end
        end
        for q = 1:K-K_FWA
            if (q~=k)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
            end
        end
        snr_num_Cell(k,n) = DS_Cell(k,n);
        snr_den_Cell(k,n) = MSI_Cell(k,n) + MUI_Cell(k,n) + HI_Cell(k,n) + noise_Cell(k,n);
        rate_dl(k+K_FWA) = rate_dl(k+K_FWA) + BW*TAU_FAC*log2(1+snr_num_Cell(k,n)/snr_den_Cell(k,n));
    end
end
end
