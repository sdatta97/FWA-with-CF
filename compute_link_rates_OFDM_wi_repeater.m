function rate_dl = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor; %scalar, or per-UE vector (per-class CSI overheads)
if isscalar(TAU_FAC)
    TAU_FAC = TAU_FAC*ones(K - params.numCPE,1);
end
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
N_CPE_FWA = size(channel_dl_FWA,4);
p_d = params.rho_tot; %total sector power. Constant-PSD convention: channels
%are normalized by FULL-band noise, so p_d*|h|^2 is the PSD ratio - exact
%on any subband; the bandwidth share enters only through the rate prelog
D_FWA = params.D_FWA;
D_cell = params.D_cell;
BETA_interUE = params.BETA_interUE;
%NCR forward-chain gain: computed AFTER attachment below, because the
%NCR forwards SUBBAND-SELECTIVELY (see there)
%Number of spatial layers the NCR forwards. 1 = Rel-18 NCR (analog
%beamforming, single amplify-and-forward chain, no digital processing:
%Carvalho et al., "Network-controlled repeater - an introduction", IEEE
%Comm. Standards Mag. 2026). Values > 1 model a beyond-Rel-18
%multi-chain "MIMO-capable" repeater forwarding the top singular beam
%pairs with the output power split equally across chains; the effective
%rank is clamped to what the arrays and the UE can support, so setting
%params.ncr_rank = params.N_UE_cell forwards as many layers as the
%cellular UE has antennas (no hardcoded layer count)
if isfield(params,'ncr_rank')
    ncr_rank = min([params.ncr_rank, N_UE, N_CPE_FWA]);
else
    ncr_rank = 1;
end
%Donor-side sector antenna of the network-controlled repeaters: amplitude
%gain of the CPE donor-locked sector panel towards sector BS m
rep_donor_amp = sqrt(params.rep_donor_gain);
%BS-to-CPE large-scale gain (linear), from the shared gain matrix; used
%only for the repeater-selection metric
BETA_FWA = params.gainOverNoise_lin(:,1:K_FWA).*D_FWA;
if params.HW_IMPAIRMENTS
    Kt = params.Kt;
    Kr = params.Kr;
    K_rep_hw = params.Kt_rep*params.Kr_rep; %repeater AF-chain quality factor
else
    Kt = 1; Kr = 1; K_rep_hw = 1;
end
%Serving sector of each UE and per-sector scheduled sets; each sector
%shares the full band equally among its UEs (OFDMA)
Serv = cell(K-K_FWA,1);
for k = 1:K-K_FWA
    Serv{k} = find(D_cell(:,k)==1);
end
U = cell(M_sectors,1);
for m = 1:M_sectors
    U{m} = find(D_cell(m,:)==1);
end
sectorLoad = sum(D_cell,2);
K_rep = params.num_repeater_per_cpe;
%EVERY served cellular UE may attach to a repeater: the per-UE benefit
%gate below keeps only genuinely helped UEs on the repeater path, so
%universal eligibility is safe by construction (strong, interference-
%limited UEs fail the gate and keep the pure direct path)
assist = false(K-K_FWA,1);
for m = 1:M_sectors
    assist(U{m}) = true;
end
%Donor sector of each CPE/repeater: the sector it is served by, which is
%where its donor panel is locked. An NCR forwards its donor gNB's signal,
%so a UE may only be assisted by a repeater whose donor sector is the
%UE's own serving sector - attaching to a repeater fed by another sector
%would inject that sector's signal as amplified interference
donor_of = zeros(K_FWA,1);
for q = 1:K_FWA
    donor_of(q) = find(D_FWA(:,q)==1,1);
end
%Repeaters (CPEs) associated with each assisted UE
Rep = cell(K-K_FWA,1);
for k = 1:K-K_FWA
    if assist(k)
        %repeater selection metric includes both repeater antenna gains:
        %donor (in rep_donor_gain) and service (in BETA_interUE)
        v = BETA_FWA(Serv{k},:)'.*params.rep_donor_gain(Serv{k},:)'.*BETA_interUE(1:K_FWA,k+K_FWA);
        %only CPEs enabled as repeaters (params.set_repeat) whose donor
        %sector matches the UE's serving sector may be selected; each
        %assisted UE associates with at most K_rep of them
        v(setdiff(1:K_FWA,params.set_repeat)) = -Inf;
        v(donor_of ~= Serv{k}) = -Inf;
        [vmax,servingCPEs] =  maxk(v,K_rep);
        Rep{k} = servingCPEs(isfinite(vmax));
    end
end
%SUBBAND-SELECTIVE FORWARDING (TR 38.867: side control switches the
%NCR-Fwd ON/OFF and steers beams PER RESOURCE): the NCR amplifies ONLY
%its attached UEs' OFDMA subbands, not the whole carrier. Consequences:
% (i) the 33 dBm output cap applies to the forwarded subbands alone -
%     with constant-PSD input, the input power on a fraction f of the
%     band is f*P_in, so the effective gain min(G_max, P_max/(f*P_in))
%     gains 10log10(1/f) dB of headroom over full-band forwarding
%     (~ +10 dB at typical sector loads);
% (ii) non-attached UEs' subbands are simply not forwarded, so only
%     attached UEs carry repeater path components - which is exactly
%     how H_rep is structured below.
n_att = zeros(K_FWA,1); %attached (repeater, UE) pairs per repeater
for k = 1:K-K_FWA
    for kk = 1:numel(Rep{k})
        n_att(Rep{k}(kk)) = n_att(Rep{k}(kk)) + 1;
    end
end
G_max_lin = 10^(0.1*params.repeat_gain);
P_max_mW = 10^(0.1*params.rep_max_pow_dBm);
rep_gain_amp = zeros(K_FWA,1);
for q = 1:K_FWA
    if n_att(q) > 0
        f_sub = n_att(q)/sectorLoad(donor_of(q)); %forwarded band fraction
        rep_gain_amp(q) = sqrt(min(G_max_lin, P_max_mW/(params.rep_P_in_mW(q)*f_sub)));
    end
end
%Repeater-path component of every sector-to-UE channel, computed once.
%The NCR applies rank-1 amplify-and-forward BEAMFORMING (TR 38.867 side
%control enables "transmissions and receptions with improved spatial
%directivity"): a donor RECEIVE beam w_r matched to its donor-sector
%link, and a service TRANSMIT beam w_t matched to the assisted user's
%link (per-user on that user's OFDMA resources). Both beams are selected
%from the ESTIMATED channels and applied to the true and estimated
%composites alike. The donor beam gives the donor sector full array gain
%while other sectors' signals land on an (almost) uncorrelated beam, so
%the repeater's input SIR improves by the array factor instead of the
%repeater re-radiating the raw interference mix:
%H_rep(m,k,:,n) = sum over UE k's repeaters of
%  (BS m -> repeater)*w_r x donor gain x repeater gain x w_t'(rep -> UE k)
H_rep = zeros(M_sectors,K-K_FWA,N_BS,N_UE);
%estimated counterpart, used only for the gNB beamforming (imperfect CSI)
H_rep_est = zeros(M_sectors,K-K_FWA,N_BS,N_UE);
%per-repeater donor receive beams, computed ONCE per attached repeater
W_r_rep = cell(K_FWA,1);
for k = 1:K-K_FWA
    for kk = 1:numel(Rep{k})
        rep_idx = Rep{k}(kk);
        if ~isempty(W_r_rep{rep_idx}), continue, end
        m_donor = donor_of(rep_idx);
        Hd_est = reshape(channel_est_dl_FWA(m_donor,rep_idx,:,:),[N_BS,N_CPE_FWA]);
        if ~isfield(params,'ncr_zf_donor') || params.ncr_zf_donor
            %ZERO-FORCING-STYLE (max-SINR / regularized ZF) donor receive
            %combining: the repeated path's SIR ceiling is set at the
            %donor input, where the co-channel sectors' transmissions
            %arrive alongside the donor's and are amplified identically.
            %The receive beams are therefore chosen as the top
            %generalized eigenvectors of (donor-subspace power,
            %interference-plus-noise covariance), suppressing the
            %interfering sites in the combining itself instead of only
            %favouring the donor. Covariances use the ESTIMATED channels
            %and the same per-sector donor-panel gains as the composite.
            Bq = eye(N_CPE_FWA); %captured noise (noise-normalized units)
            for mm = setdiff(1:M_sectors,m_donor)
                Hi = reshape(channel_est_dl_FWA(mm,rep_idx,:,:),[N_BS,N_CPE_FWA]);
                Bq = Bq + (p_d/N_BS)*params.rep_donor_gain(mm,rep_idx)*(Hi'*Hi);
            end
            [Vg,Dg] = eig(Hd_est'*Hd_est, Bq, 'vector');
            [~,ordg] = sort(real(Dg),'descend');
            Vr = Vg(:,ordg);
            Vr = Vr./vecnorm(Vr); %unit-norm receive beams per chain
        else
            %matched-only donor combining (ablation): dominant donor
            %directions, no interference awareness
            [~,~,Vr] = svd(Hd_est,'econ');
        end
        W_r_rep{rep_idx} = Vr(:,1:min(ncr_rank,size(Vr,2)));
    end
end
for k = 1:K-K_FWA
    for kk = 1:numel(Rep{k})
        rep_idx = Rep{k}(kk);
        Vr = W_r_rep{rep_idx};
        Hs_est = reshape(channel_interFWA_est(rep_idx,k+K_FWA,:,1:N_UE),[N_CPE_FWA,N_UE]);
        [~,~,Vt] = svd(Hs_est.','econ');
        %forwarding chains: strongest donor direction pairs with strongest
        %service direction; the amplifier output (33 dBm cap) splits
        %equally across chains, hence the 1/sqrt(ncr_rank) per-chain
        %amplitude
        for chain = 1:ncr_rank
            w_r = Vr(:,chain); %donor receive beam of this chain
            w_t = Vt(:,chain); %service transmit beam of this chain
            amp_chain = rep_gain_amp(rep_idx)/sqrt(ncr_rank);
            for m = 1:M_sectors
                g_vec = reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,rep_idx)*amp_chain;
                g_vec_est = reshape(channel_est_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,rep_idx)*amp_chain;
                for n = 1:N_UE
                    s_n = reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[1,N_CPE_FWA])*w_t;
                    s_n_est = reshape(channel_interFWA_est(rep_idx,k+K_FWA,:,n),[1,N_CPE_FWA])*w_t;
                    H_rep(m,k,:,n) = reshape(H_rep(m,k,:,n),[N_BS,1]) + g_vec*s_n;
                    H_rep_est(m,k,:,n) = reshape(H_rep_est(m,k,:,n),[N_BS,1]) + g_vec_est*s_n_est;
                end
            end
        end
    end
end
H_comp = channel_dl + H_rep; %repeater-inclusive composite channels (true)
H_comp_est = channel_est_dl + H_rep_est; %estimated composites, for beamforming only
%Matched-filter beams each sector applies towards its own scheduled UEs,
%normalized columns of the ESTIMATED composite channels; physical
%reception uses the true ones. The same matrices provide the victim's own
%stream beams and the beamformed inter-sector interference, averaged over
%the interferer's scheduled set (round-robin OFDMA)
W = cell(M_sectors,1);
W_dir = cell(M_sectors,1); %direct-channel beams, for the NCR-OFF variant
for m = 1:M_sectors
    W{m} = beamMatrix(H_comp_est,m,U{m},N_BS,N_UE);
    W_dir{m} = beamMatrix(channel_est_dl,m,U{m},N_BS,N_UE);
end

%% Computing rates: MIMO log-det per UE
%Shannon capacity of UE k's N_UE x N_UE effective channel A (receive
%antennas x streams) under colored interference-plus-noise:
%  rate_k = share*prelog*log2 det(I + Kt*Kr*A' R_k^{-1} A),
%R_k = beamformed inter-sector interference as a spatial covariance
%(averaged over the interferer's scheduled beams, round-robin OFDMA)
%+ diagonal EVM/repeater distortion + noise. The NCR-ON/OFF benefit
%gate compares the two log-det rates.
noise_dl = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    m = Serv{k};
    share = BW/sectorLoad(m);
    W_own = W{m}; U_own = U{m};
    W_dir_own = W_dir{m};
    int_sectors = setdiff(1:M_sectors,m); %reuse 1: every other sector interferes
    %NCR-OFF (pure direct path) rate, evaluated for every UE: for UEs
    %without an attached repeater this IS their rate; for attached UEs it
    %is the side-control alternative. With the NCR not forwarding on this
    %UE's resources there is no re-amplified interference either.
    Wk_dir = zeros(N_BS,N_UE);
    for n = 1:N_UE
        Wk_dir(:,n) = getBeam(W_dir_own,U_own,k,n,N_UE);
    end
    H_dir = reshape(channel_dl(m,k,:,1:N_UE),[N_BS,N_UE]); %columns = receive antennas
    %Per-sector radiated-power constraint: p_d TOTAL on the subband, split
    %evenly over the N_UE unit-norm stream beams (matching the /(|U_m|*N_UE)
    %divisor in the interference averages)
    A_off = sqrt(p_d/N_UE)*(H_dir'*Wk_dir);
    R_mci_off = zeros(N_UE);
    for mm = int_sectors
        if ~isempty(W{mm})
            H_int = reshape(channel_dl(mm,k,:,1:N_UE),[N_BS,N_UE]);
            B = sqrt(p_d)*(H_int'*W{mm});
            R_mci_off = R_mci_off + (B*B')/size(W{mm},2);
        end
    end
    R_off = diag(noise_dl(k,:)) + Kr*Kt*R_mci_off ...
          + diag((1-Kr*Kt)*(sum(abs(A_off).^2,2) + real(diag(R_mci_off))));
    rate_off = share*TAU_FAC(k)*logdet2(Kr*Kt*(A_off'*(R_off\A_off)));
    rate_dl(k) = rate_off;
    if ~isempty(Rep{k})
        %NCR-ON (repeater composite) rate. Benefit-gated side control
        %(TR 38.867): the gNB keeps the better of ON/OFF per UE, since an
        %AF repeater substitutes a rank-1, interference-inheriting path
        %for the direct MIMO link and only pays off in coverage holes.
        %params.ncr_benefit_gate = 0 disables the gate (ablation): every
        %attached UE is forced onto the repeater composite.
        %Cross-UE coupling through the interferers' beam matrices is
        %neglected in the OFF evaluation (measured < 1% of a UE's rate).
        %The useful signal keeps only sqrt(K_rep_hw) of the repeater
        %branch; the rest becomes repeater hardware distortion (diagonal,
        %uncorrelated-distortion convention, matching the scalar model).
        Wk_comp = zeros(N_BS,N_UE); %the sector's beams on the composite channel
        for n = 1:N_UE
            Wk_comp(:,n) = getBeam(W_own,U_own,k,n,N_UE);
        end
        H_rep_k = reshape(H_rep(m,k,:,1:N_UE),[N_BS,N_UE]);
        A_sig = sqrt(p_d/N_UE)*((H_dir + sqrt(K_rep_hw)*H_rep_k)'*Wk_comp);
        A_rep = sqrt(p_d/N_UE)*(H_rep_k'*Wk_comp);
        %inter-sector interference with the interferers' actual
        %beamforming; the victim's interference channel includes its
        %repeater re-amplifying the interfering sector's signal
        R_mci_on = zeros(N_UE);
        for mm = int_sectors
            if ~isempty(W{mm})
                H_int = reshape(H_comp(mm,k,:,1:N_UE),[N_BS,N_UE]);
                B = sqrt(p_d)*(H_int'*W{mm});
                R_mci_on = R_mci_on + (B*B')/size(W{mm},2);
            end
        end
        R_on = diag(noise_dl(k,:)) + Kr*Kt*R_mci_on ...
             + diag((1-Kr*Kt)*(sum(abs(A_sig).^2,2) + real(diag(R_mci_on))) ...
                  + (1-K_rep_hw)*sum(abs(A_rep).^2,2));
        rate_on = share*TAU_FAC(k)*logdet2(Kr*Kt*(A_sig'*(R_on\A_sig)));
        gate_on = ~isfield(params,'ncr_benefit_gate') || params.ncr_benefit_gate;
        if rate_on > rate_off || ~ gate_on
            rate_dl(k) = rate_on;
        end
    end
end
end

function W = beamMatrix(H_comp,m,U,N_BS,N_UE)
%Columns are the unit-norm matched-filter beams sector m applies towards
%its scheduled UEs (one per UE antenna stream), from the composite
%estimated channels
if isempty(U)
    W = zeros(N_BS,0);
    return
end
W = zeros(N_BS,numel(U)*N_UE);
c = 0;
for q = U(:)'
    for nn = 1:N_UE
        c = c + 1;
        h = reshape(H_comp(m,q,:,nn),N_BS,1);
        W(:,c) = h./norm(h);
    end
end
end

function w = getBeam(W,U,q,nn,N_UE)
%Beam a sector applies towards antenna nn of its scheduled UE q: a column
%of its beam matrix
w = W(:,(find(U==q,1)-1)*N_UE+nn);
end
