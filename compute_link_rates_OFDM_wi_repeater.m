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
rep_gain = params.repeat_gain;
%Donor-side sector antenna of the network-controlled repeaters: amplitude
%gain of CPE q's best-pointing sector panel towards sector BS m
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
%Repeaters (CPEs) associated with each UE
Rep = cell(K-K_FWA,1);
for k = 1:K-K_FWA
    %repeater selection metric includes both repeater antenna gains:
    %donor (in rep_donor_gain) and service (in BETA_interUE)
    v = BETA_FWA(Serv{k},:)'.*params.rep_donor_gain(Serv{k},:)'.*BETA_interUE(1:K_FWA,k+K_FWA);
    %only CPEs enabled as repeaters (params.set_repeat) may be selected;
    %each UE associates with at most K_rep of them
    v(setdiff(1:K_FWA,params.set_repeat)) = -Inf;
    [vmax,servingCPEs] =  maxk(v,K_rep);
    Rep{k} = servingCPEs(isfinite(vmax));
end
%Repeater-path component of every sector-to-UE channel, computed once:
%H_rep(m,k,:,n) = sum over UE k's repeaters of (BS m -> repeater) x donor
%gain x repeater gain x (repeater -> UE k). The composite channel used
%for beams, DS/MSI, and interference is channel_dl + H_rep.
H_rep = zeros(M_sectors,K-K_FWA,N_BS,N_UE);
%estimated counterpart, built from the ESTIMATED BS-to-repeater and
%repeater-to-UE channels; used only for beamforming (imperfect CSI)
H_rep_est = zeros(M_sectors,K-K_FWA,N_BS,N_UE);
for k = 1:K-K_FWA
    for kk = 1:numel(Rep{k})
        rep_idx = Rep{k}(kk);
        for m = 1:M_sectors
            G = reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain;
            G_est = reshape(channel_est_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain;
            for n = 1:N_UE
                H_rep(m,k,:,n) = reshape(H_rep(m,k,:,n),[N_BS,1]) + G*reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[N_CPE_FWA,1]);
                H_rep_est(m,k,:,n) = reshape(H_rep_est(m,k,:,n),[N_BS,1]) + G_est*reshape(channel_interFWA_est(rep_idx,k+K_FWA,:,n),[N_CPE_FWA,1]);
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
for m = 1:M_sectors
    W{m} = beamMatrix(H_comp_est,m,U{m},N_BS,N_UE);
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
    int_sectors = setdiff(1:M_sectors,m); %reuse 1: every other sector interferes
    for n = 1:N_UE
        %receive channel of UE k split into direct and repeater-path
        %components: the useful signal keeps only sqrt(K_rep_hw) of the
        %repeater branch, the rest becomes repeater hardware distortion
        eff_dir = reshape(channel_dl(m,k,:,n),N_BS,1);
        eff_rep = reshape(H_rep(m,k,:,n),N_BS,1);
        eff_channel = eff_dir + eff_rep;            %composite (true)
        eff_sig = eff_dir + sqrt(K_rep_hw)*eff_rep; %coherent part after repeater HI
        w_n = getBeam(W_own,U_own,k,n,N_UE);
        ds_base = p_d*abs(eff_sig'*w_n)^2;
        DS_dl(k,n) = Kr*Kt*ds_base;
        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base ...
                   + (1-K_rep_hw)*p_d*abs(eff_rep'*w_n)^2;
        for nn = 1:N_UE
            if (nn~=n)
                nn_eff_channel = reshape(H_comp(m,k,:,nn),N_BS,1);
                if (norm(nn_eff_channel,'fro') < norm(eff_channel,'fro'))
                    w_nn = getBeam(W_own,U_own,k,nn,N_UE);
                    msi_base = p_d*abs(eff_sig'*w_nn)^2;
                    MSI_dl(k,n) = Kr*Kt*msi_base;
                    HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base ...
                               + (1-K_rep_hw)*p_d*abs(eff_rep'*w_nn)^2;
                end
            end
        end
        %inter-sector (co-channel) interference with the interferers'
        %actual beamforming towards their own scheduled UEs; the victim's
        %interference channel includes its repeaters re-amplifying the
        %interfering sector's signal (donor gain towards that sector)
        mci_base = 0;
        for mm = int_sectors
            if ~isempty(W{mm})
                h_int = reshape(H_comp(mm,k,:,n),N_BS,1);
                mci_base = mci_base + p_d*mean(abs(h_int'*W{mm}).^2);
            end
        end
        MCI_dl(k,n) = Kr*Kt*mci_base;
        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base;
        rate_dl(k) = rate_dl(k) + share*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
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
