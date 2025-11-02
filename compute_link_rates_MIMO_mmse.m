function rate_dl = compute_link_rates_MIMO_mmse(params,channel_dl, channel_dl_sc, channel_est_dl, channel_est_dl_sc, channel_dl_FWA, channel_dl_FWA_sc, channel_est_dl_FWA, channel_est_dl_FWA_sc)
M = params.numGNB;
S = params.numSC;
K_FWA = params.numCPE;
K = params.numUE + params.numCPE;
BW = params.Band;
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_SC = size(channel_dl_sc,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
p_d_sc = params.rho_tot_sc; % 1*K;
p_fac = params.p_fac;
D = params.D;
BETA = params.BETA;
D_sc = params.D_sc;
BETA_sc = params.BETA_sc;
beta_uc = zeros(size(BETA));
beta_uc_sc = zeros(size(BETA_sc));

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv_sc = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ_sc = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingBSs = find(D(:,k)==1);
    NoservingBSs = find(D(:,k)==0);
    servingSCs = find(D_sc(:,k)==1);
    NoservingSCs = find(D_sc(:,k)==0);    
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
    Serv_sc{k} = servingSCs;
    NoServ_sc{k} = NoservingSCs;  

    La(k) = length(servingBSs);
    beta_uc(:,k) = BETA(:,k).*D(:,k);
    beta_uc_sc(:,k) = BETA_sc(:,k).*D_sc(:,k);
end

%% initialization of c
D_FWA_FWA = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
D_FWA_Cell = zeros(K_FWA,K-K_FWA,N_CPE_FWA,N_UE);
D_Cell_FWA = zeros(K-K_FWA,K_FWA,N_UE,N_CPE_FWA);
D_Cell_Cell = zeros(K-K_FWA,K-K_FWA,N_UE,N_UE);
dl_mmse_precoder_FWA = zeros(size(channel_est_dl_FWA));
dl_mmse_precoder = zeros(size(channel_est_dl));
scaling_LP_mmse = zeros(M,K);
dl_mmse_precoder_FWA_sc = zeros(size(channel_est_dl_FWA_sc));
dl_mmse_precoder_sc = zeros(size(channel_est_dl_sc));
scaling_LP_mmse_sc = zeros(S,K);
for m = 1:M
    for k = 1:K_FWA
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(N_BS);
        for q = 1:K_FWA
            if ismember(m,Serv{q})                 
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])';
            end
        end
        for q = 1:K-K_FWA
            if ismember(m,Serv{q+K_FWA})
                inv_matrix = inv_matrix + p_d*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])';
            end
        end
        dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA]) + p_d*inv_matrix\(reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]));
        if ismember(m,Serv{k})
            scaling_LP_mmse(m,k) = scaling_LP_mmse(m,k) + norm(dl_mmse_precoder_FWA(m,k,:,:),'fro')^2;
        end
    end
    for k = 1:K-K_FWA
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(N_BS);
        for q = 1:K_FWA
            if ismember(m,Serv{q}) 
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*reshape(channel_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])';
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
for s = 1:S
    for k = 1:K_FWA
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(N_SC);
        for q = 1:K_FWA
            if ismember(s,Serv_sc{q})                 
                inv_matrix = inv_matrix + p_d_sc*reshape(channel_dl_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA])*reshape(channel_dl_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA])';
            end
        end
        for q = 1:K-K_FWA
            if ismember(s,Serv_sc{q+K_FWA})
                inv_matrix = inv_matrix + p_d_sc*reshape(channel_dl_sc(s,q,:,:),[N_SC,N_UE])*reshape(channel_dl_sc(s,q,:,:),[N_SC,N_UE])';
            end
        end
        dl_mmse_precoder_FWA_sc(s,k,:,:) = reshape(dl_mmse_precoder_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA]) + p_d_sc*inv_matrix\(reshape(channel_dl_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA]));
        if ismember(s,Serv_sc{k})
            scaling_LP_mmse_sc(s,k) = scaling_LP_mmse_sc(s,k) + norm(dl_mmse_precoder_FWA_sc(s,k,:,:),'fro')^2;
        end
    end
    for k = 1:K-K_FWA
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(N_SC);
        for q = 1:K_FWA
            if ismember(s,Serv_sc{q}) 
                inv_matrix = inv_matrix + p_d_sc*reshape(channel_dl_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA])*reshape(channel_dl_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA])';
            end
        end
        for q = 1:K-K_FWA
            if ismember(s,Serv_sc{q+K_FWA})
                inv_matrix = inv_matrix +  p_d_sc*reshape(channel_dl_sc(s,q,:,:),[N_SC,N_UE])*reshape(channel_dl_sc(s,q,:,:),[N_SC,N_UE])';
            end
        end
        dl_mmse_precoder_sc(s,k,:,:) = reshape(dl_mmse_precoder_sc(s,k,:,:),[N_SC,N_UE]) + p_d_sc*inv_matrix\(reshape(channel_dl_sc(s,k,:,:),[N_SC,N_UE]));
        if ismember(s,Serv_sc{k+K_FWA})
            scaling_LP_mmse_sc(s,k+K_FWA) = scaling_LP_mmse_sc(s,k+K_FWA) + norm(dl_mmse_precoder_sc(s,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M
    for k = 1:K_FWA
        if ismember(m,Serv{k})
            dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])./sqrt(scaling_LP_mmse(m,k));
        end
    end
    for k = 1:K-K_FWA
        if ismember(m,Serv{k+K_FWA})
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE])./sqrt(scaling_LP_mmse(m,k+K_FWA));
        end
    end
end
for s = 1:S
    for k = 1:K_FWA
        if ismember(s,Serv_sc{k})
            dl_mmse_precoder_FWA_sc(s,k,:,:) = reshape(dl_mmse_precoder_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA])./sqrt(scaling_LP_mmse_sc(s,k));
        end
    end
    for k = 1:K-K_FWA
        if ismember(s,Serv_sc{k+K_FWA})
            dl_mmse_precoder_sc(s,k,:,:) = reshape(dl_mmse_precoder_sc(s,k,:,:),[N_SC,N_UE])./sqrt(scaling_LP_mmse_sc(s,k+K_FWA));
        end
    end
end
eta_eq = zeros(M,K);
for m = 1:M
    term = 0;
    for k = 1:K
        if ismember(m,Serv{k})
            if (k<=K_FWA)
                term = term + p_fac*trace(reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])');
            else
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])*reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])');
            end
        end
    end
    if (term > 0)
        eta_eq(m,:) = (1/term)*D(m,:);
        eta_eq(m,1:K_FWA) = p_fac*eta_eq(m,1:K_FWA);
    end
end
eta_eq_sc = zeros(S,K);
for s = 1:S
    term = 0;
    for k = 1:K
        if ismember(s,Serv_sc{k})
            if (k<=K_FWA)
                term = term + p_fac*trace(reshape(dl_mmse_precoder_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA])');
            else
                term = term + trace(reshape(dl_mmse_precoder_sc(m,k-K_FWA,:,:),[N_SC,N_UE])*reshape(dl_mmse_precoder_sc(s,k-K_FWA,:,:),[N_SC,N_UE])');
            end
        end
    end
    if (term > 0)
        eta_eq_sc(m,:) = (1/term)*D_sc(m,:);
        eta_eq_sc(m,1:K_FWA) = p_fac*eta_eq_sc(m,1:K_FWA);
    end
end
for k = 1:K_FWA
    for q = 1:K_FWA
        for m = 1:M
            if ismember(m,Serv{q})
                D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(eta_eq(m,q))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
            end
        end
        for s = 1:S
            if ismember(s,Serv_sc{q})
                D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(eta_eq_sc(s,q))*reshape(channel_dl_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA])'*reshape(dl_mmse_precoder_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA]);
            end
        end
    end
    for q = 1:K-K_FWA
        for m = 1:M
            if ismember(m,Serv{q+K_FWA})
                D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(eta_eq(m,q+K_FWA))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
            end
        end
        for s = 1:S
            if ismember(s,Serv_sc{q+K_FWA})
                D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(eta_eq_sc(s,q+K_FWA))*reshape(channel_dl_FWA_sc(s,k,:,:),[N_SC,N_CPE_FWA])'*reshape(dl_mmse_precoder_sc(s,q,:,:),[N_SC,N_UE]);
            end
        end
    end
end
for k = 1:K-K_FWA
    for q = 1:K_FWA
        for m = 1:M
            if ismember(m,Serv{q})
                D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
            end
        end
        for s = 1:S
            if ismember(s,Serv_sc{q})
                D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d_sc*eta_eq_sc(s,q))*reshape(channel_dl_sc(s,k,:,:),[N_SC,N_UE])'*reshape(dl_mmse_precoder_FWA_sc(s,q,:,:),[N_SC,N_CPE_FWA]);
            end
        end
    end
    for q = 1:K-K_FWA
        for m = 1:M
            if ismember(m,Serv{q+K_FWA})
                D_Cell_Cell(k,q,:,:) = reshape(D_Cell_Cell(k,q,:,:),[N_UE,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
            end
        end
        for s = 1:S
            if ismember(s,Serv_sc{q+K_FWA})
                D_Cell_Cell(k,q,:,:) = reshape(D_Cell_Cell(k,q,:,:),[N_UE,N_UE]) + sqrt(p_d_sc*eta_eq_sc(s,q+K_FWA))*reshape(channel_dl_sc(s,k,:,:),[N_SC,N_UE])'*reshape(dl_mmse_precoder_sc(s,q,:,:),[N_SC,N_UE]);
            end
        end
    end
end
DS_FWA = zeros(K_FWA,N_CPE_FWA);
MSI_FWA = zeros(K_FWA,N_CPE_FWA);
MUI_FWA = zeros(K_FWA,N_CPE_FWA);
DS_Cell = zeros(K-K_FWA,N_UE);
MSI_Cell = zeros(K-K_FWA,N_UE);
MUI_Cell = zeros(K-K_FWA,N_UE);

noise_FWA = abs(sqrt(0.5)*(randn(K_FWA,N_CPE_FWA) + 1j*randn(K_FWA,N_CPE_FWA))).^2;
noise_Cell = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
snr_num_FWA = zeros(K_FWA,N_CPE_FWA);
snr_den_FWA = zeros(K_FWA,N_CPE_FWA);
snr_num_Cell = zeros(K-K_FWA,N_UE);
snr_den_Cell = zeros(K-K_FWA,N_UE);
rate_dl = zeros(K,1);
for k = 1:K_FWA
    for n = 1:N_CPE_FWA
        DS_FWA(k,n) = (abs(D_FWA_FWA(k,k,n,n)))^2;
        for nn = 1:N_CPE_FWA
            if (abs(D_FWA_FWA(k,k,nn,nn))<abs(D_FWA_FWA(k,k,n,n)))
                MSI_FWA(k,n) = MSI_FWA(k,n) + (abs(D_FWA_FWA(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_FWA
            if (q~=k)
              MUI_FWA(k,n) = MUI_FWA(k,n) + norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
            end
        end
        for q = 1:K-K_FWA
           MUI_FWA(k,n) = MUI_FWA(k,n) + norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
        end
        snr_num_FWA(k,n) = DS_FWA(k,n);
        snr_den_FWA(k,n) = MSI_FWA(k,n) + MUI_FWA(k,n) + noise_FWA(k,n);
        rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_FWA(k,n)/snr_den_FWA(k,n));
    end
end
for k = 1:K-K_FWA
    for n = 1:N_UE
        DS_Cell(k,n) = (abs(D_Cell_Cell(k,k,n,n)))^2;
        for nn = 1:N_UE
            if (abs(D_Cell_Cell(k,k,nn,nn))<abs(D_Cell_Cell(k,k,n,n)))
                MSI_Cell(k,n) = MSI_Cell(k,n) + (abs(D_Cell_Cell(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_FWA
            if (q~=k)
              MUI_Cell(k,n) = MUI_Cell(k,n) + norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
            end
        end
        for q = 1:K-K_FWA
           MUI_Cell(k,n) = MUI_Cell(k,n) + norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
        end
        snr_num_Cell(k,n) = DS_Cell(k,n);
        snr_den_Cell(k,n) = MSI_Cell(k,n) + MUI_Cell(k,n) + noise_Cell(k,n);
        rate_dl(k+K_FWA) = rate_dl(k+K_FWA) + BW*TAU_FAC*log2(1+snr_num_Cell(k,n)/snr_den_Cell(k,n));
    end
end
end