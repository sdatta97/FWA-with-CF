function rate_dl = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
N_CPE_FWA = size(channel_dl_FWA,4);
p_d = params.rho_tot; % 1*K;
D_FWA = params.D_FWA;
D_cell = params.D_cell;
BETA_interUE = params.BETA_interUE;
%Per-repeater effective AMPLITUDE gain of the NCR forward chain: fixed
%amplification gain capped by the maximum output power at the repeater's
%operating point (TR 38.867: 90 dB / 33 dBm), computed in the main script
rep_gain_amp = params.repeat_gain_amp;
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
%Repeater assistance is restricted to the WEAKEST cellular UEs of each
%sector (the cell-edge users): the rep_assist_frac fraction with the
%lowest serving-link large-scale gain. Strong UEs keep the pure direct
%path, so repeater-injected interference cannot degrade the healthy
%majority (and starve the cell band adaptation at its low percentile).
assist = false(K-K_FWA,1);
gain_ue = params.gainOverNoise_lin(:,K_FWA+1:end);
for m = 1:M_sectors
    served = U{m};
    [~,ord] = sort(gain_ue(m,served),'ascend');
    n_assist = ceil(params.rep_assist_frac*numel(served));
    assist(served(ord(1:n_assist))) = true;
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
for k = 1:K-K_FWA
    for kk = 1:numel(Rep{k})
        rep_idx = Rep{k}(kk);
        m_donor = donor_of(rep_idx);
        Hd_est = reshape(channel_est_dl_FWA(m_donor,rep_idx,:,:),[N_BS,N_CPE_FWA]);
        [~,~,Vr] = svd(Hd_est,'econ');
        w_r = Vr(:,1);   %donor receive beam (dominant donor-link direction)
        Hs_est = reshape(channel_interFWA_est(rep_idx,k+K_FWA,:,1:N_UE),[N_CPE_FWA,N_UE]);
        [~,~,Vt] = svd(Hs_est.','econ');
        w_t = Vt(:,1);   %service transmit beam towards assisted UE k
        for m = 1:M_sectors
            g_vec = reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,rep_idx)*rep_gain_amp(rep_idx);
            g_vec_est = reshape(channel_est_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*w_r*rep_donor_amp(m,rep_idx)*rep_gain_amp(rep_idx);
            for n = 1:N_UE
                s_n = reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[1,N_CPE_FWA])*w_t;
                s_n_est = reshape(channel_interFWA_est(rep_idx,k+K_FWA,:,n),[1,N_CPE_FWA])*w_t;
                H_rep(m,k,:,n) = reshape(H_rep(m,k,:,n),[N_BS,1]) + g_vec*s_n;
                H_rep_est(m,k,:,n) = reshape(H_rep_est(m,k,:,n),[N_BS,1]) + g_vec_est*s_n_est;
            end
        end
    end
end
H_comp = channel_dl + H_rep; %repeater-inclusive composite channels (true)
H_comp_est = channel_est_dl + H_rep_est; %estimated composites, for beamforming only
%Matched-filter beams each sector applies towards its own scheduled UEs,
%normalized columns of the ESTIMATED composite channels; physical
%reception uses the true ones. The same matrices provide the victim's own
%DS/MSI beams and the beamformed inter-sector interference, averaged over
%the interferer's scheduled set (round-robin OFDMA)
W = cell(M_sectors,1);
W_dir = cell(M_sectors,1); %direct-channel beams, for the NCR-OFF variant
for m = 1:M_sectors
    W{m} = beamMatrix(H_comp_est,m,U{m},N_BS,N_UE);
    W_dir{m} = beamMatrix(channel_est_dl,m,U{m},N_BS,N_UE);
end

%% Computing rates
DS_dl = zeros(K-K_FWA,N_UE);
MSI_dl = zeros(K-K_FWA,N_UE);
MCI_dl = zeros(K-K_FWA,N_UE);
HI_dl = zeros(K-K_FWA,N_UE);
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
    rate_off = 0;
    for n = 1:N_UE
        h_n = reshape(channel_dl(m,k,:,n),N_BS,1);
        w_n = getBeam(W_dir_own,U_own,k,n,N_UE);
        ds_base = p_d*abs(h_n'*w_n)^2;
        DS_off = Kr*Kt*ds_base;
        HI_off = (1-Kr*Kt)*ds_base;
        MSI_off = 0;
        for nn = 1:N_UE
            if (nn~=n)
                h_nn = reshape(channel_dl(m,k,:,nn),N_BS,1);
                if (h_nn'*h_nn < h_n'*h_n)
                    msi_base = p_d*abs(h_n'*getBeam(W_dir_own,U_own,k,nn,N_UE))^2;
                    MSI_off = Kr*Kt*msi_base;
                    HI_off = HI_off + (1-Kr*Kt)*msi_base;
                end
            end
        end
        mci_base = 0;
        for mm = int_sectors
            if ~isempty(W{mm})
                h_int = reshape(channel_dl(mm,k,:,n),N_BS,1);
                mci_base = mci_base + p_d*mean(abs(h_int'*W{mm}).^2);
            end
        end
        MCI_off = Kr*Kt*mci_base;
        HI_off = HI_off + (1-Kr*Kt)*mci_base;
        DS_dl(k,n) = DS_off; MSI_dl(k,n) = MSI_off;
        MCI_dl(k,n) = MCI_off; HI_dl(k,n) = HI_off;
        rate_off = rate_off + share*TAU_FAC*log2(1+DS_off/(MSI_off+MCI_off+HI_off+noise_dl(k,n)));
    end
    rate_dl(k) = rate_off;
    if ~isempty(Rep{k})
        %NCR-ON (repeater composite) rate. Benefit-gated side control
        %(TR 38.867): the gNB keeps the better of ON/OFF per UE, since an
        %AF repeater substitutes a rank-1, interference-inheriting path
        %for the direct MIMO link and only pays off in coverage holes.
        %Cross-UE coupling through the interferers' beam matrices is
        %neglected in the OFF evaluation (measured < 1% of a UE's rate).
        rate_on = 0;
        DS_on = zeros(1,N_UE); MSI_on = zeros(1,N_UE);
        MCI_on = zeros(1,N_UE); HI_on = zeros(1,N_UE);
        for n = 1:N_UE
            %receive channel split into direct and repeater-path parts:
            %the useful signal keeps only sqrt(K_rep_hw) of the repeater
            %branch, the rest becomes repeater hardware distortion
            eff_dir = reshape(channel_dl(m,k,:,n),N_BS,1);
            eff_rep = reshape(H_rep(m,k,:,n),N_BS,1);
            eff_channel = eff_dir + eff_rep;            %composite (true)
            eff_sig = eff_dir + sqrt(K_rep_hw)*eff_rep; %coherent part after repeater HI
            w_n = getBeam(W_own,U_own,k,n,N_UE);
            ds_base = p_d*abs(eff_sig'*w_n)^2;
            DS_on(n) = Kr*Kt*ds_base;
            HI_on(n) = HI_on(n) + (1-Kr*Kt)*ds_base ...
                     + (1-K_rep_hw)*p_d*abs(eff_rep'*w_n)^2;
            for nn = 1:N_UE
                if (nn~=n)
                    nn_eff_channel = reshape(H_comp(m,k,:,nn),N_BS,1);
                    if (norm(nn_eff_channel,'fro') < norm(eff_channel,'fro'))
                        w_nn = getBeam(W_own,U_own,k,nn,N_UE);
                        msi_base = p_d*abs(eff_sig'*w_nn)^2;
                        MSI_on(n) = Kr*Kt*msi_base;
                        HI_on(n) = HI_on(n) + (1-Kr*Kt)*msi_base ...
                                 + (1-K_rep_hw)*p_d*abs(eff_rep'*w_nn)^2;
                    end
                end
            end
            %inter-sector interference with the interferers' actual
            %beamforming; the victim's interference channel includes its
            %repeater re-amplifying the interfering sector's signal
            mci_base = 0;
            for mm = int_sectors
                if ~isempty(W{mm})
                    h_int = reshape(H_comp(mm,k,:,n),N_BS,1);
                    mci_base = mci_base + p_d*mean(abs(h_int'*W{mm}).^2);
                end
            end
            MCI_on(n) = Kr*Kt*mci_base;
            HI_on(n) = HI_on(n) + (1-Kr*Kt)*mci_base;
            rate_on = rate_on + share*TAU_FAC*log2(1+DS_on(n)/(MSI_on(n)+MCI_on(n)+HI_on(n)+noise_dl(k,n)));
        end
        if rate_on > rate_off
            rate_dl(k) = rate_on;
            DS_dl(k,:) = DS_on; MSI_dl(k,:) = MSI_on;
            MCI_dl(k,:) = MCI_on; HI_dl(k,:) = HI_on;
        end
    end
end
end

function W = beamMatrix(H_comp,m,U,N_BS,N_UE)
%Columns are the unit-norm matched-filter beams sector m applies towards
%its scheduled UEs (one per UE antenna stream), from the composite channels
if isempty(U)
    W = zeros(N_BS,0);
    return
end
W = zeros(N_BS,numel(U)*N_UE);
c = 0;
for q = U(:)'
    for nn = 1:N_UE
        h = reshape(H_comp(m,q,:,nn),N_BS,1);
        c = c + 1;
        W(:,c) = h./norm(h);
    end
end
end

function w = getBeam(W,U,q,nn,N_UE)
%Beam a sector applies towards antenna nn of its scheduled UE q: a column
%of its beam matrix
w = W(:,(find(U==q,1)-1)*N_UE+nn);
end
