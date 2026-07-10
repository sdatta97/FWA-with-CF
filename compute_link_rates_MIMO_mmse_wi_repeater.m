function rate_dl = compute_link_rates_MIMO_mmse_wi_repeater(params,channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D_FWA = params.D_FWA;
%BS-to-CPE large-scale gain (linear), from the shared gain matrix; used
%only for the repeater-selection metric
BETA = params.gainOverNoise_lin(:,1:K_FWA);
BETA_interUE = params.BETA_interUE;
rep_gain = params.repeat_gain;
%Donor-side sector antenna of the network-controlled repeaters: amplitude
%gain of CPE q's best-pointing sector panel towards sector BS m
rep_donor_amp = sqrt(params.rep_donor_gain);
BETA = BETA.*D_FWA;
if params.HW_IMPAIRMENTS
    Kt      = params.Kt;
    Kr_FWA  = params.Kr;
    Kr_cell = params.Kr;
    K_rep_hw = params.Kt_rep*params.Kr_rep; %repeater AF-chain quality factor
else
    Kt = 1; Kr_FWA = 1; Kr_cell = 1; K_rep_hw = 1;
end
%Serving sector(s) of each user
Serv = cell(K,1);
for k = 1:K
    Serv{k} = find(D_FWA(:,k)==1);
end

K_rep = params.num_repeater_per_cpe;
%Repeaters (CPEs) associated with each served CPE, in the same format as
%the OFDM rate functions: one selection per user, capped at K_rep
Rep = cell(K_FWA,1);
for k = 1:K_FWA
    if ~ismember(k,params.set_repeat)
        %repeater selection metric includes both repeater antenna gains:
        %donor (in rep_donor_gain) and service (in BETA_interUE),
        %aggregated over the CPE's serving sector(s)
        v = sum(BETA(Serv{k},1:K_FWA).*params.rep_donor_gain(Serv{k},:),1)'.*BETA_interUE(1:K_FWA,k);
        %only CPEs enabled as repeaters (params.set_repeat) may be
        %selected; each served CPE associates with at most K_rep of them
        v(setdiff(1:K_FWA,params.set_repeat)) = -Inf;
        [vmax,servingCPEs] =  maxk(v,K_rep);
        Rep{k} = sort(servingCPEs(isfinite(vmax)));
    end
end
%Repeater-path component of every sector-to-CPE channel, computed once
%(analogous to the composite tensor of the OFDM rate functions):
%H_FWA_rep(m,k,:,:) = sum over CPE k's repeaters of (BS m -> repeater) x
%donor gain x repeater gain x (repeater -> CPE k). The composite channel
%used for the MMSE filters, effective gains, and interference is
%channel_dl_FWA + H_FWA_rep.
H_FWA_rep = zeros(M_sectors,K_FWA,N_BS,N_CPE_FWA);
%estimated counterpart, built from the ESTIMATED BS-to-repeater and
%repeater-to-CPE channels; used only for beamforming (imperfect CSI)
H_FWA_rep_est = zeros(M_sectors,K_FWA,N_BS,N_CPE_FWA);
for k = 1:K_FWA
    for kk = 1:numel(Rep{k})
        rep_idx = Rep{k}(kk);
        for m = 1:M_sectors
            H_FWA_rep(m,k,:,:) = reshape(H_FWA_rep(m,k,:,:),[N_BS,N_CPE_FWA]) + reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,k,:,:),[N_CPE_FWA,N_CPE_FWA]);
            H_FWA_rep_est(m,k,:,:) = reshape(H_FWA_rep_est(m,k,:,:),[N_BS,N_CPE_FWA]) + reshape(channel_est_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA_est(rep_idx,k,:,:),[N_CPE_FWA,N_CPE_FWA]);
        end
    end
end
%% initialization of c
D_FWA_FWA = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
D_FWA_Cell = zeros(K_FWA,K-K_FWA,N_CPE_FWA,N_UE);
%repeater-path-only components of the effective gains, used to compute
%the distortion added by the repeater hardware (fraction 1-K_rep_hw of the
%power arriving through a repeater)
D_FWA_FWA_rep = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
D_FWA_Cell_rep = zeros(K_FWA,K-K_FWA,N_CPE_FWA,N_UE);
D_Cell_FWA = zeros(K-K_FWA,K_FWA,N_UE,N_CPE_FWA);
D_Cell_Cell = zeros(K-K_FWA,K-K_FWA,N_UE,N_UE);
dl_mmse_precoder_FWA = zeros(size(channel_est_dl_FWA));
dl_mmse_precoder = zeros(size(channel_est_dl));
scaling_LP_mmse = zeros(M_sectors,K);
%MMSE precoding, structured like the OFDM beam matrices: the MMSE filter
%basis inv_matrix is identical for every user a sector precodes, so it
%is built once per sector (on the repeater-inclusive composite channels)
%and applied to each served user's channel. Beamforming operates on the
%ESTIMATED composites (imperfect CSI when params.IMPERFECT_CSI is set);
%physical reception below uses the true ones
for m = 1:M_sectors
    inv_matrix = eye(N_BS);
    for q = 1:K_FWA
        if ismember(m,Serv{q}) && ~ismember(q,params.set_repeat)
            eff_channel = reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA]) + reshape(H_FWA_rep_est(m,q,:,:),[N_BS,N_CPE_FWA]);
            inv_matrix = inv_matrix + p_d*(eff_channel*eff_channel');
        end
    end
    for q = 1:K-K_FWA
        if ismember(m,Serv{q+K_FWA})
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])';
        end
    end
    for k = 1:K_FWA
        if ~ismember(k,params.set_repeat)
            eff_channel = reshape(channel_est_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]) + reshape(H_FWA_rep_est(m,k,:,:),[N_BS,N_CPE_FWA]);
            dl_mmse_precoder_FWA(m,k,:,:) = p_d*inv_matrix\(eff_channel);
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
        if ismember(m,Serv{k}) && ~ismember(k,params.set_repeat)
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
                if ~ismember(k,params.set_repeat)
                    term = term + trace(reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])');
                end
            else
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])*reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])');
            end
        end
    end
    if (term > 0)
        setA = 1:K;
        eta_eq(m,setA(~ismember(setA,params.set_repeat))) = (1/term)*D_FWA(m,setA(~ismember(setA,params.set_repeat)));
    end
end
%Effective gains: the channels already carry the sector antenna pattern
%of BS m towards each user (via the pattern-inclusive large-scale gains
%inside R_gNB) and the repeater donor-panel gains, so both the desired
%links and every inter-sector interference term below are sectored
for k = 1:K_FWA
    if ~ismember(k,params.set_repeat)
        %Effective receive channel of CPE k from each sector BS m, split
        %into the direct component and the repeater-path component
        %(tensor lookups). The useful signal carries only sqrt(K_rep_hw)
        %of the repeater branch; the rest of the repeater-path power
        %becomes hardware distortion
        for q = 1:K_FWA
            if ~ismember(q,params.set_repeat)
                for m = 1:M_sectors
                    if ismember(m,Serv{q})
                        eff_dir = reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]);
                        eff_rep = reshape(H_FWA_rep(m,k,:,:),[N_BS,N_CPE_FWA]);
                        w_q = reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                        D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*(eff_dir + sqrt(K_rep_hw)*eff_rep)'*w_q;
                        D_FWA_FWA_rep(k,q,:,:) = reshape(D_FWA_FWA_rep(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*eff_rep'*w_q;
                    end
                end
            end
        end
        for q = 1:K-K_FWA
            for m = 1:M_sectors
                if ismember(m,Serv{q+K_FWA})
                    eff_dir = reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]);
                    eff_rep = reshape(H_FWA_rep(m,k,:,:),[N_BS,N_CPE_FWA]);
                    w_q = reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
                    D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*(eff_dir + sqrt(K_rep_hw)*eff_rep)'*w_q;
                    D_FWA_Cell_rep(k,q,:,:) = reshape(D_FWA_Cell_rep(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*eff_rep'*w_q;
                end
            end
        end
    end
end
for k = 1:K-K_FWA
    for q = 1:K_FWA
        if ~ismember(q,params.set_repeat)
            for m = 1:M_sectors
                if ismember(m,Serv{q})
                    D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                end
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

noise_FWA = abs(sqrt(0.5)*(randn(K_FWA,N_CPE_FWA) + 1j*randn(K_FWA,N_CPE_FWA))).^2;
noise_Cell = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
snr_num_FWA = zeros(K_FWA,N_CPE_FWA);
snr_den_FWA = zeros(K_FWA,N_CPE_FWA);
snr_num_Cell = zeros(K-K_FWA,N_UE);
snr_den_Cell = zeros(K-K_FWA,N_UE);
rate_dl = zeros(K,1);
for k = 1:K_FWA
    if ~ismember(k,params.set_repeat)
        for n = 1:N_CPE_FWA
            DS_FWA(k,n) = Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,n)))^2;
            HI_FWA(k,n) = HI_FWA(k,n) + (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,n)))^2;
            for nn = 1:N_CPE_FWA
                if (abs(D_FWA_FWA(k,k,nn,nn))<abs(D_FWA_FWA(k,k,n,n)))
                    MSI_FWA(k,n) = MSI_FWA(k,n) + Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                    HI_FWA(k,n) = HI_FWA(k,n) + (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                end
            end
            %repeater hardware distortion from the own-signal streams:
            %fraction (1-K_rep_hw) of the own-stream power that reaches CPE
            %k through a repeater path
            HI_FWA(k,n) = HI_FWA(k,n) + (1-K_rep_hw)*norm(reshape(D_FWA_FWA_rep(k,k,n,:),[1,N_CPE_FWA]))^2;
            for q = 1:K_FWA
                if (q~=k) && ~ismember(q,params.set_repeat)
                  MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                  HI_FWA(k,n) = HI_FWA(k,n) + (1-Kr_FWA*Kt)*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2 ...
                              + (1-K_rep_hw)*norm(reshape(D_FWA_FWA_rep(k,q,n,:),[1,N_CPE_FWA]))^2;
                end
            end
            for q = 1:K-K_FWA
               MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
               HI_FWA(k,n) = HI_FWA(k,n) + (1-Kr_FWA*Kt)*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2 ...
                           + (1-K_rep_hw)*norm(reshape(D_FWA_Cell_rep(k,q,n,:),[1,N_UE]))^2;
            end
            snr_num_FWA(k,n) = DS_FWA(k,n);
            snr_den_FWA(k,n) = MSI_FWA(k,n) + MUI_FWA(k,n) + HI_FWA(k,n) + noise_FWA(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_FWA(k,n)/snr_den_FWA(k,n));
        end
    end
end
for k = 1:K-K_FWA
    for n = 1:N_UE
        DS_Cell(k,n) = Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,n)))^2;
        HI_Cell(k,n) = HI_Cell(k,n) + (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,n)))^2;
        for nn = 1:N_UE
            if (abs(D_Cell_Cell(k,k,nn,nn))<abs(D_Cell_Cell(k,k,n,n)))
                MSI_Cell(k,n) = MSI_Cell(k,n) + Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,nn)))^2;
                HI_Cell(k,n) = HI_Cell(k,n) + (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_FWA
            if ~ismember(q,params.set_repeat)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                HI_Cell(k,n) = HI_Cell(k,n) + (1-Kr_cell*Kt)*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
            end
        end
        for q = 1:K-K_FWA
            if (q~=k)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
                HI_Cell(k,n) = HI_Cell(k,n) + (1-Kr_cell*Kt)*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
            end
        end
        snr_num_Cell(k,n) = DS_Cell(k,n);
        snr_den_Cell(k,n) = MSI_Cell(k,n) + MUI_Cell(k,n) + HI_Cell(k,n) + noise_Cell(k,n);
        rate_dl(k+K_FWA) = rate_dl(k+K_FWA) + BW*TAU_FAC*log2(1+snr_num_Cell(k,n)/snr_den_Cell(k,n));
    end
end
end
