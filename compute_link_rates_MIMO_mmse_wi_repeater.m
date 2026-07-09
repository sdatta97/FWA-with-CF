function rate_dl = compute_link_rates_MIMO_mmse_wi_repeater(params,channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est)
M = params.numGNB;
K_FWA = params.numCPE;
K = M*params.numUE + params.numCPE;
BW = params.Band;
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D_FWA = params.D_FWA;
BETA = params.BETA;
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
%Prepare cell to store the AP indices serving a specific UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specific UE
NoServ = cell(K,1);
%Prepare cell to store the AP indices serving a specific UE
%Prepare cell to store the AP indices not serving a specific UE
%Construct the above array and cells
for k = 1:K
    servingBSs = find(D_FWA(:,k)==1);
    NoservingBSs = find(D_FWA(:,k)==0);
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
end

K_rep = params.num_repeater_per_cpe;
%Prepare cell to store the AP indices serving a specific UE
Rep = cell(K_FWA,1);
%Prepare cell to store the AP indices not serving a specific UE
NoRep = cell(K_FWA,1);
%Prepare cell to store the AP indices serving a specific UE
%Prepare cell to store the AP indices not serving a specific UE
%Construct the above array and cells
for k = 1:K_FWA
    if ~ismember(k,params.set_repeat)
        for m = 1:M
            if ismember(m,Serv{k})
                %repeater selection metric includes both repeater antenna
                %gains: donor (in rep_donor_gain) and service (in BETA_interUE)
                v = BETA(m,1:K_FWA)'.*params.rep_donor_gain(m,1:K_FWA)'.*BETA_interUE(1:K_FWA,k);
                v(setdiff(1:K_FWA,params.set_repeat)) = -Inf;
                [vmax,servingCPEs] =  maxk(v,K_rep);
                servingCPEs = servingCPEs(isfinite(vmax)); %guard: pool may be smaller than K_rep
                NoservingCPEs = setdiff(params.set_repeat,servingCPEs);
                Rep{k} = union(Rep{k},servingCPEs);
                NoRep{k} = union(NoRep{k}, NoservingCPEs); 
            end
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
scaling_LP_mmse = zeros(M,K);
for m = 1:M
    for k = 1:K_FWA
        if ~ismember(k,params.set_repeat)
            inv_matrix = eye(N_BS);
            for q = 1:K_FWA
                if ismember(m,Serv{q}) && ~ismember(q,params.set_repeat) 
                    eff_channel = reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                    for qq = 1:numel(Rep{q})
                        rep_idx = Rep{q}(qq);
                        eff_channel = eff_channel + reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,q,:,:),[N_CPE_FWA,N_CPE_FWA]);
                    end
                    inv_matrix = inv_matrix + p_d*(eff_channel*eff_channel');
                end
            end
            for q = 1:K-K_FWA
                if ismember(m,Serv{q+K_FWA})
                    inv_matrix = inv_matrix + p_d*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])';
                end
            end
            eff_channel = reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]);
            for kk = 1:numel(Rep{k})
                rep_idx = Rep{k}(kk);
                eff_channel = eff_channel + reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,k,:,:),[N_CPE_FWA,N_CPE_FWA]);
            end
            dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA]) + p_d*inv_matrix\(eff_channel);
            if ismember(m,Serv{k})
                scaling_LP_mmse(m,k) = scaling_LP_mmse(m,k) + norm(dl_mmse_precoder_FWA(m,k,:,:),'fro')^2;
            end
        end
    end
    for k = 1:K-K_FWA
        inv_matrix = eye(N_BS);
        for q = 1:K_FWA
            if ismember(m,Serv{q}) && ~ismember(q,params.set_repeat)
                eff_channel = reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                for qq = 1:numel(Rep{q})
                    rep_idx = Rep{q}(qq);
                    eff_channel = eff_channel + reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,q,:,:),[N_CPE_FWA,N_CPE_FWA]);
                end
                inv_matrix = inv_matrix + p_d*(eff_channel*eff_channel');
            end
        end
        for q = 1:K-K_FWA
            if ismember(m,Serv{q+K_FWA})
                inv_matrix = inv_matrix +  p_d*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])';
            end
        end
        dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE]) + p_d*inv_matrix\(reshape(channel_dl(m,k,:,:),[N_BS,N_UE]));
        if ismember(m,Serv{k+K_FWA})
            scaling_LP_mmse(m,k+K_FWA) = scaling_LP_mmse(m,k+K_FWA) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M
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
eta_eq = zeros(M,K);
for m = 1:M
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
for k = 1:K_FWA
    if ~ismember(k,params.set_repeat)
        %Effective receive channel of CPE k from each sector BS m, split
        %into the direct component and the repeater-path component. The
        %useful signal carries only sqrt(K_rep_hw) of the repeater branch;
        %the rest of the repeater-path power becomes hardware distortion
        eff_dir_all = cell(M,1);
        eff_rep_all = cell(M,1);
        for m = 1:M
            eff_dir_all{m} = reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]);
            eff_rep_all{m} = zeros(N_BS,N_CPE_FWA);
            for qq = 1:numel(Rep{k})
                rep_idx = Rep{k}(qq);
                eff_rep_all{m} = eff_rep_all{m} + reshape(channel_dl_FWA(m,rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(m,rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,k,:,:),[N_CPE_FWA,N_CPE_FWA]);
            end
        end
        for q = 1:K_FWA
            if ~ismember(q,params.set_repeat)
                for m = 1:M
                    if ismember(m,Serv{q})
                        w_q = reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                        D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*(eff_dir_all{m} + sqrt(K_rep_hw)*eff_rep_all{m})'*w_q;
                        D_FWA_FWA_rep(k,q,:,:) = reshape(D_FWA_FWA_rep(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*eff_rep_all{m}'*w_q;
                    end
                end
            end
        end
        for q = 1:K-K_FWA
            for m = 1:M
                if ismember(m,Serv{q+K_FWA})
                    w_q = reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
                    D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*(eff_dir_all{m} + sqrt(K_rep_hw)*eff_rep_all{m})'*w_q;
                    D_FWA_Cell_rep(k,q,:,:) = reshape(D_FWA_Cell_rep(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*eff_rep_all{m}'*w_q;
                end
            end
        end
    end
end
for k = 1:K-K_FWA
    for q = 1:K_FWA
        if ~ismember(q,params.set_repeat) 
            for m = 1:M
                if ismember(m,Serv{q})
                    D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                end
            end
        end
    end
    for q = 1:K-K_FWA
        for m = 1:M
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