function rate_dl = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW)
M = params.numGNB;
K_mmW = params.numUE;
K = params.numUE + params.numUE_sub6;
BW = params.Band;
TAU_FAC = params.preLogFactor;
Ntx = size(channel_dl,3);
N_UE_mmW = size(channel_dl_mmW,4);
N_UE_sub6 = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D = params.D;
BETA = params.BETA;
beta_uc = zeros(size(BETA));

%Prepare array to store the number of APs serving a specficic UE
La = zeros(K,1);
%Prepare cell to store the AP indices serving a specficic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specficic UE
NoServ = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingAPs = find(D(:,k)==1);
    NoservingAPs = find(D(:,k)==0);
    
    Serv{k} = servingAPs;
    NoServ{k} = NoservingAPs;
    
    La(k) = length(servingAPs);
    beta_uc(:,k) = BETA(:,k).*D(:,k);
end

%% initialization of c
D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_mmW,N_UE_mmW);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_mmW,N_UE_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_sub6,N_UE_mmW);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_sub6,N_UE_sub6);
dl_mmse_precoder_mmW = zeros(size(channel_est_dl_mmW));
dl_mmse_precoder = zeros(size(channel_est_dl));
scaling_LP_mmse = zeros(M,K);
for m = 1:M
    for k = 1:K_mmW
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(Ntx);
        for q = 1:K_mmW
            if ismember(m,Serv{q})                 
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])';
            end
        end
        for q = 1:K-K_mmW
            if ismember(m,Serv{q+K_mmW})
                inv_matrix = inv_matrix + p_d*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])';
            end
        end
        dl_mmse_precoder_mmW(m,k,:,:) = reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW]) + p_d*inv_matrix\(reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW]));
        if ismember(m,Serv{k})
            scaling_LP_mmse(m,k) = scaling_LP_mmse(m,k) + norm(dl_mmse_precoder_mmW(m,k,:,:),'fro')^2;
        end
    end
    for k = 1:K-K_mmW
        % inv_matrix = noiseVariance*eye(Ntx);
        inv_matrix = eye(Ntx);
        for q = 1:K_mmW
            if ismember(m,Serv{q}) 
                inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[Ntx,N_UE_mmW])';
            end
        end
        for q = 1:K-K_mmW
            if ismember(m,Serv{q+K_mmW})
                inv_matrix = inv_matrix +  p_d*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[Ntx,N_UE_sub6])';
            end
        end
        dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[Ntx,N_UE_sub6]) + p_d*inv_matrix\(reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6]));
        if ismember(m,Serv{k+K_mmW})
            scaling_LP_mmse(m,k+K_mmW) = scaling_LP_mmse(m,k+K_mmW) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M
    for k = 1:K_mmW
        if ismember(m,Serv{k})
            dl_mmse_precoder_mmW(m,k,:,:) = reshape(dl_mmse_precoder_mmW(m,k,:,:),[Ntx,N_UE_mmW])./sqrt(scaling_LP_mmse(m,k));
        end
    end
    for k = 1:K-K_mmW
        if ismember(m,Serv{k+K_mmW})
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[Ntx,N_UE_sub6])./sqrt(scaling_LP_mmse(m,k+K_mmW));
        end
    end
end
eta_eq = zeros(M,K);
N_AP = params.num_antennas_per_gNB;
if (K_mmW == 0)
    for m = 1:M
        term = 0;
        for k = 1+K_mmW:K
            if ismember(m,Serv{k})
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[N_AP,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[N_AP,N_UE_sub6])');
            end
        end
        if (term > 0)
            eta_eq(m,:) = (1/term)*D(m,:);
        end
    end
else
    for m = 1:M
        term = 0;
        for k = 1:K
            if ismember(m,Serv{k})
                if (k<=K_mmW)
                    term = term + trace(reshape(dl_mmse_precoder_mmW(m,k,:,:),[N_AP,N_UE_mmW])*reshape(dl_mmse_precoder_mmW(m,k,:,:),[N_AP,N_UE_mmW])');
                else
                    term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[N_AP,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:),[N_AP,N_UE_sub6])');
                end
            end
        end
        if (term > 0)
            eta_eq(m,:) = (1/term)*D(m,:);
        end
    end
end
for k = 1:K_mmW
    for q = 1:K_mmW
        for m = 1:M
            if ismember(m,Serv{q})
                D_mmW_mmW(k,q,:,:) = reshape(D_mmW_mmW(k,q,:,:),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
            end
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            if ismember(m,Serv{q+K_mmW})
                D_mmW_sub6(k,q,:,:) = reshape(D_mmW_sub6(k,q,:,:),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(channel_dl_mmW(m,k,:,:),[Ntx,N_UE_mmW])'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
            end
        end
    end
end
for k = 1:K-K_mmW
    for q = 1:K_mmW
        for m = 1:M
            if ismember(m,Serv{q})
                D_sub6_mmW(k,q,:,:) = reshape(D_sub6_mmW(k,q,:,:),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder_mmW(m,q,:,:),[Ntx,N_UE_mmW]);
            end
        end
    end
    for q = 1:K-K_mmW
        for m = 1:M
            if ismember(m,Serv{q+K_mmW})
                D_sub6_sub6(k,q,:,:) = reshape(D_sub6_sub6(k,q,:,:),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW))*reshape(channel_dl(m,k,:,:),[Ntx,N_UE_sub6])'*reshape(dl_mmse_precoder(m,q,:,:),[Ntx,N_UE_sub6]);
            end
        end
    end
end
DS_mmW = zeros(K_mmW,N_UE_mmW);
MSI_mmW = zeros(K_mmW,N_UE_mmW);
MUI_mmW = zeros(K_mmW,N_UE_mmW);
DS_sub6 = zeros(K-K_mmW,N_UE_sub6);
MSI_sub6 = zeros(K-K_mmW,N_UE_sub6);
MUI_sub6 = zeros(K-K_mmW,N_UE_sub6);

noise_mmW = abs(sqrt(0.5)*(randn(K_mmW,N_UE_mmW) + 1j*randn(K_mmW,N_UE_mmW))).^2;
noise_sub6 = abs(sqrt(0.5)*(randn(K-K_mmW,N_UE_sub6) + 1j*randn(K-K_mmW,N_UE_sub6))).^2;
snr_num_mmW = zeros(K_mmW,N_UE_mmW);
snr_den_mmW = zeros(K_mmW,N_UE_mmW);
snr_num_sub6 = zeros(K-K_mmW,N_UE_sub6);
snr_den_sub6 = zeros(K-K_mmW,N_UE_sub6);
rate_dl = zeros(K,1);
for k = 1:K_mmW
    for n = 1:N_UE_mmW
        DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n)))^2;
        for nn = 1:N_UE_mmW
            if (abs(D_mmW_mmW(k,k,nn,nn))<abs(D_mmW_mmW(k,k,n,n)))
                MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_mmW
            if (q~=k)
              MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
            end
        end
        for q = 1:K-K_mmW
           MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
        end
        snr_num_mmW(k,n) = DS_mmW(k,n);
        snr_den_mmW(k,n) = MSI_mmW(k,n) + MUI_mmW(k,n) + noise_mmW(k,n);
        rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
    end
end
for k = 1:K-K_mmW
    for n = 1:N_UE_sub6
        DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n)))^2;
        for nn = 1:N_UE_sub6
            if (abs(D_sub6_sub6(k,k,nn,nn))<abs(D_sub6_sub6(k,k,n,n)))
                MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_mmW
            if (q~=k)
              MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:),[1,N_UE_mmW]))^2;
            end
        end
        for q = 1:K-K_mmW
           MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:),[1,N_UE_sub6]))^2;
        end
        snr_num_sub6(k,n) = DS_sub6(k,n);
        snr_den_sub6(k,n) = MSI_sub6(k,n) + MUI_sub6(k,n) + noise_sub6(k,n);
        rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + BW*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
    end
end
end