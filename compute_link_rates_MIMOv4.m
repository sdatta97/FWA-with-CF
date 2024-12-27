function rate_dl = compute_link_rates_MIMOv4(params,channel_dl, channel_est_dl,channel_dl_mmW, channel_est_dl_mmW,sub6ConnectionState)
M = size(channel_dl,1);
K_mmW = size(sub6ConnectionState,1);
K = K_mmW + size(channel_dl,2);
BW = params.Band;
scs = params.scs_sub6;
num_sc_sub6 = params.num_sc_sub6;
K_sc = K/num_sc_sub6;
TAU_FAC = params.preLogFactor;
N_AP = params.num_antennas_per_gNB;
N_UE_FWA = params.N_UE_FWA;
N_UE_cell = params.N_UE_cell;
p_d = params.rho_tot; % 1*K;
p_fac = params.p_fac;
p_fac_rearrange = params.p_fac_rearrange;
D = params.D;
% ue_rearranged = params.ue_rearranged;
user_sc_alloc = params.user_sc_alloc;
% perm_vec  = repmat(randperm(tau_p),1,2);
% phi_index = perm_vec(1:K);
% for k = 1:K
%     PHI(:,k) = PHI1(:,phi_index(k));
% end
% gamma_num = zeros(M,K);
% gamma_den = zeros(M,K);
% Gamma = zeros(M,K);
% for m = 1:M
%     for k = 1:K
%         gamma_num(m,k) = tau_p*p_p*(BETA(m,k)^2);
%         gamma_den_temp = zeros(1,K);
%         for j = 1:K
%             gamma_den_temp(j) = BETA(m,j)*(abs(PHI(:,j)'*PHI(:,k))^2);
%         end
% 
%         gamma_den(m,k) = tau_p*p_p*sum(gamma_den_temp)+1;
% 
%         Gamma(m,k) = gamma_num(m,k)/gamma_den(m,k);
%     end
% end
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

% dl_mmse_precoder_mmW = zeros([size(channel_est_dl_mmW),num_sc_sub6]);
% dl_mmse_precoder = zeros([size(channel_est_dl),num_sc_sub6]);
dl_mmse_precoder_mmW = zeros([M,K_mmW,N_AP,N_UE_FWA,num_sc_sub6]);
dl_mmse_precoder = zeros([M,K-K_mmW,N_AP,N_UE_cell,num_sc_sub6]);
scaling_LP_mmse = zeros(M,K,num_sc_sub6);
for m = 1:M
    for n = 1:num_sc_sub6
        for k = 1:K_mmW
            if (sub6ConnectionState(k) == 1) && (user_sc_alloc(k,n) == 1)
                % inv_matrix = noiseVariance*eye(N_AP);
                inv_matrix = eye(N_AP);
                for q = 1:K_mmW
                    if ismember(m,Serv{q}) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_mmW])';
                        inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_FWA])*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_FWA])';
                    end
                end
                for q = 1:K-K_mmW
                    if ismember(m,Serv{q+K_mmW}) && (user_sc_alloc(q+K_mmW,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_sub6])';
                        inv_matrix = inv_matrix + p_d*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_cell])*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_cell])';
                    end
                end
                dl_mmse_precoder_mmW(m,k,:,:,n) = reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_FWA]) + p_d*inv_matrix\(reshape(channel_dl_mmW(m,k,:,:),[N_AP,N_UE_FWA]));
                if ismember(m,Serv{k})
                    % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder_mmW(m,k,:,:),'fro')^2;
                    scaling_LP_mmse(m,k,n) = scaling_LP_mmse(m,k,n) + norm(dl_mmse_precoder_mmW(m,k,:,:,n),'fro')^2;
                end
            end
        end
        for k = 1:K-K_mmW
            if (user_sc_alloc(k+K_mmW,n) == 1)
                % inv_matrix = noiseVariance*eye(N_AP);
                inv_matrix = eye(N_AP);
                for q = 1:K_mmW
                    if ismember(m,Serv{q}) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n) == 1)
                        % inv_matrix = inv_matrix + p_d*noiseVariance*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_mmW])*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_mmW])';
                        inv_matrix = inv_matrix + p_d*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_FWA])*reshape(channel_dl_mmW(m,q,:,:),[N_AP,N_UE_FWA])';
                    end
                end
                for q = 1:K-K_mmW
                    if ismember(m,Serv{q+K_mmW}) && (user_sc_alloc(q+K_mmW,n) == 1)
                        % inv_matrix = inv_matrix +  p_d*noiseVariance*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_sub6])*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_sub6])';
                        inv_matrix = inv_matrix +  p_d*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_cell])*reshape(channel_dl(m,q,:,:),[N_AP,N_UE_cell])';
                    end
                end
                dl_mmse_precoder(m,k,:,:,n) = reshape(dl_mmse_precoder(m,k,:,:,n),[N_AP,N_UE_cell]) + p_d*inv_matrix\(reshape(channel_dl(m,k,:,:),[N_AP,N_UE_cell]));
                if ismember(m,Serv{k+K_mmW})
                    % scaling_LP_mmse(m,servedUEs) = scaling_LP_mmse(m,servedUEs) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
                    scaling_LP_mmse(m,k+K_mmW,n) = scaling_LP_mmse(m,k+K_mmW,n) + norm(dl_mmse_precoder(m,k,:,:,n),'fro')^2;
                end
            end
        end
    end
end
% dl_mmse_precoder_mmW = conj(channel_est_dl_mmW);
% dl_mmse_precoder = conj(channel_est_dl);
for m = 1:M
    for n = 1:num_sc_sub6
        for k = 1:K_mmW
            if ismember(m,Serv{k}) && (user_sc_alloc(k,n) == 1) && (sub6ConnectionState(k)==1)
                dl_mmse_precoder_mmW(m,k,:,:,n) = reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_FWA])./sqrt(scaling_LP_mmse(m,k,n));
            end
        end
        for k = 1:K-K_mmW
            if ismember(m,Serv{k+K_mmW}) && (user_sc_alloc(k+K_mmW,n) == 1)
                dl_mmse_precoder(m,k,:,:,n) = reshape(dl_mmse_precoder(m,k,:,:,n),[N_AP,N_UE_cell])./sqrt(scaling_LP_mmse(m,k+K_mmW,n));
            end
        end
    end
end

%% initialization of c
eta_eq = zeros(M,K,num_sc_sub6);
if (K_mmW == 0) || all(sub6ConnectionState == 0)
    for m = 1:M
        for n = 1:num_sc_sub6
            term = 0;
            for k = 1+K_mmW:K
                if ismember(m,Serv{k}) && (user_sc_alloc(k,n) == 1)
                        %term = (N_AP*N_UE_sub6*num_sc_sub6*beta_uc(m,:)*user_sc_alloc(:,n));
                        % term = num_sc_sub6*N_AP*(N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                      
                        % term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_sub6])');
                        term = term + num_sc_sub6*trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_cell])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_cell])');
                end
            end
            if (term > 0)
                eta_eq(m,:,n) = (1/term)*(D(m,:)'.*(user_sc_alloc(:,n)==1));
            end
        end
    end
else
    % ues_not_rearranged = setdiff((1+K_mmW):K,ue_rearranged);
    for m = 1:M
        for n = 1:num_sc_sub6
            term = 0;
            for k = 1:K
                if ismember(m,Serv{k}) && (user_sc_alloc(k,n) == 1)
                    if (k<=K_mmW) && (sub6ConnectionState(k) == 1) 
    %                         term = (N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                        % term = num_sc_sub6*N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                                              
                        % term = term + p_fac*trace(reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_mmW])*reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_mmW])');
                        term = term + num_sc_sub6*p_fac*trace(reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_FWA])*reshape(dl_mmse_precoder_mmW(m,k,:,:,n),[N_AP,N_UE_FWA])');
                    elseif (k>K_mmW)
%                     eta_eq(m,k) = 1./(N_AP*(N_UE_mmW*p_fac*beta_uc(m,1:K_mmW)+N_UE_sub6*sum(beta_uc(m,2:K))));
%                         term = (N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*beta_uc(m,(1+K_mmW):K)*user_sc_alloc((1+K_mmW):K,n)));
                        % term = num_sc_sub6*N_AP*(N_UE_mmW*num_sc_sub6*p_fac*((beta_uc(m,1:K_mmW).*sub6ConnectionState)*user_sc_alloc(1:K_mmW,n))+N_UE_sub6*(p_fac_rearrange*beta_uc(m,ue_rearranged)*user_sc_alloc(ue_rearranged,n)+beta_uc(m,ues_not_rearranged)*user_sc_alloc(ues_not_rearranged,n)));                        
                        % term = term + trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_sub6])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_sub6])');
                        term = term + num_sc_sub6*trace(reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_cell])*reshape(dl_mmse_precoder(m,k-K_mmW,:,:,n),[N_AP,N_UE_cell])');
                    end
                end
            end
            if (term > 0)
                eta_eq(m,:,n) = (1/term)*(D(m,:)'.*(user_sc_alloc(:,n)==1));
                eta_eq(m,1:K_mmW,n) = p_fac*(eta_eq(m,1:K_mmW,n).*(sub6ConnectionState==1)');
            end
        end
    end
end
D_mmW_mmW = zeros(K_mmW,K_mmW,N_UE_FWA,N_UE_FWA,num_sc_sub6);
D_mmW_sub6 = zeros(K_mmW,K-K_mmW,N_UE_FWA,N_UE_cell,num_sc_sub6);
D_sub6_mmW = zeros(K-K_mmW,K_mmW,N_UE_cell,N_UE_FWA,num_sc_sub6);
D_sub6_sub6 = zeros(K-K_mmW,K-K_mmW,N_UE_cell,N_UE_cell,num_sc_sub6);
for n_idx = 1:num_sc_sub6
    for k = 1:K_mmW
        for q = 1:K_mmW
            for m = 1:M
                if ismember(m,Serv{q}) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n_idx) == 1)
                    % D_mmW_mmW(k,q,:,:,n_idx) = reshape(D_mmW_mmW(k,q,:,:,n_idx),[N_UE_mmW,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(channel_dl_mmW(m,k,:,:),[N_AP,N_UE_mmW]).'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[N_AP,N_UE_mmW]);
                    D_mmW_mmW(k,q,:,:,n_idx) = reshape(D_mmW_mmW(k,q,:,:,n_idx),[N_UE_FWA,N_UE_FWA]) + sqrt(eta_eq(m,q,n_idx))*reshape(channel_dl_mmW(m,k,:,:),[N_AP,N_UE_FWA])'*reshape(dl_mmse_precoder_mmW(m,q,:,:,n_idx),[N_AP,N_UE_FWA]);
                end
            end
        end
        for q = 1:K-K_mmW
            for m = 1:M
                if ismember(m,Serv{q+K_mmW}) && (user_sc_alloc(q+K_mmW,n_idx) == 1)
                    % D_mmW_sub6(k,q,:,:,n_idx) = reshape(D_mmW_sub6(k,q,:,:,n_idx),[N_UE_mmW,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(channel_dl_mmW(m,k,:,:),[N_AP,N_UE_mmW]).'*reshape(conj(channel_est_dl(m,q,:,:)),[N_AP,N_UE_sub6]);
                    D_mmW_sub6(k,q,:,:,n_idx) = reshape(D_mmW_sub6(k,q,:,:,n_idx),[N_UE_FWA,N_UE_cell]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(channel_dl_mmW(m,k,:,:),[N_AP,N_UE_FWA])'*reshape(dl_mmse_precoder(m,q,:,:,n_idx),[N_AP,N_UE_cell]);
                end
            end
        end
    end
    for k = 1:K-K_mmW
        for q = 1:K_mmW
            for m = 1:M
                if ismember(m,Serv{q}) && (sub6ConnectionState(q) == 1) && (user_sc_alloc(q,n_idx) == 1)
                    % D_sub6_mmW(k,q,:,:,n_idx) = reshape(D_sub6_mmW(k,q,:,:,n_idx),[N_UE_sub6,N_UE_mmW]) + sqrt(eta_eq(m,q,n_idx))*reshape(channel_dl(m,k,:,:),[N_AP,N_UE_sub6]).'*reshape(conj(channel_est_dl_mmW(m,q,:,:)),[N_AP,N_UE_mmW]);
                    D_sub6_mmW(k,q,:,:,n_idx) = reshape(D_sub6_mmW(k,q,:,:,n_idx),[N_UE_cell,N_UE_FWA]) + sqrt(eta_eq(m,q,n_idx))*reshape(channel_dl(m,k,:,:),[N_AP,N_UE_cell])'*reshape(dl_mmse_precoder_mmW(m,q,:,:,n_idx),[N_AP,N_UE_FWA]);
                end
            end
        end
        for q = 1:K-K_mmW
            for m = 1:M
                if ismember(m,Serv{q+K_mmW}) && (user_sc_alloc(q+K_mmW,n_idx) == 1)
                    % D_sub6_sub6(k,q,:,:,n_idx) = reshape(D_sub6_sub6(k,q,:,:,n_idx),[N_UE_sub6,N_UE_sub6]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(channel_dl(m,k,:,:),[N_AP,N_UE_sub6]).'*reshape(conj(channel_est_dl(m,q,:,:)),[N_AP,N_UE_sub6]);
                    D_sub6_sub6(k,q,:,:,n_idx) = reshape(D_sub6_sub6(k,q,:,:,n_idx),[N_UE_cell,N_UE_cell]) + sqrt(eta_eq(m,q+K_mmW,n_idx))*reshape(channel_dl(m,k,:,:),[N_AP,N_UE_cell])'*reshape(dl_mmse_precoder(m,q,:,:,n_idx),[N_AP,N_UE_cell]);
                end
            end
        end
    end
end
rate_dl = zeros(K,1);
for n_idx = 1:num_sc_sub6
    DS_mmW = zeros(K_mmW,N_UE_FWA);
    MSI_mmW = zeros(K_mmW,N_UE_FWA);
    MUI_mmW = zeros(K_mmW,N_UE_FWA);
    DS_sub6 = zeros(K-K_mmW,N_UE_cell);
    MSI_sub6 = zeros(K-K_mmW,N_UE_cell);
    MUI_sub6 = zeros(K-K_mmW,N_UE_cell);
    
    noise_mmW = abs(sqrt(0.5)*(randn(K_mmW,N_UE_FWA) + 1j*randn(K_mmW,N_UE_FWA))).^2;
    noise_sub6 = abs(sqrt(0.5)*(randn(K-K_mmW,N_UE_cell) + 1j*randn(K-K_mmW,N_UE_cell))).^2;
    snr_num_mmW = zeros(K_mmW,N_UE_FWA);
    snr_den_mmW = zeros(K_mmW,N_UE_FWA);
    snr_num_sub6 = zeros(K-K_mmW,N_UE_cell);
    snr_den_sub6 = zeros(K-K_mmW,N_UE_cell);
    for k = 1:K_mmW
    %     if (sub6ConnectionState(k)==1 || k==ue_idx)
        if (sub6ConnectionState(k)==1) && (user_sc_alloc(k,n_idx) == 1)
            for n = 1:N_UE_FWA              
    %             DS_mmW(k,n) = p_d*norm(reshape(D_mmW_mmW(k,k,n,:),[1,N_UE_mmW]))^2;
                DS_mmW(k,n) = p_d*(abs(D_mmW_mmW(k,k,n,n,n_idx)))^2;
                for nn = 1:N_UE_FWA
    %                 if (nn~=n)
                    if (abs(D_mmW_mmW(k,k,nn,nn,n_idx))<abs(D_mmW_mmW(k,k,n,n,n_idx)))
    %                     MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,k,nn,:),[1,N_UE_mmW]))^2;
                        MSI_mmW(k,n) = MSI_mmW(k,n) + p_d*(abs(D_mmW_mmW(k,k,n,nn,n_idx)))^2;
                    end
                end
                for q = 1:K_mmW
                    if (q~=k && sub6ConnectionState(q)==1)
                      MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_mmW(k,q,n,:,n_idx),[1,N_UE_FWA]))^2;
                    end
                end
                for q = 1:K-K_mmW
                   MUI_mmW(k,n) = MUI_mmW(k,n) + p_d*norm(reshape(D_mmW_sub6(k,q,n,:,n_idx),[1,N_UE_cell]))^2;
                end
                snr_num_mmW(k,n) = DS_mmW(k,n);
                snr_den_mmW(k,n) = MSI_mmW(k,n) + MUI_mmW(k,n) + noise_mmW(k,n);
                % rate_dl(k) = rate_dl(k) + scs(1)*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
                rate_dl(k) = rate_dl(k) + scs(n_idx)*TAU_FAC*log2(1+snr_num_mmW(k,n)/snr_den_mmW(k,n));
            end
        end
    end
    for k = 1:K-K_mmW
        for n = 1:N_UE_cell
           if (user_sc_alloc(k+K_mmW,n_idx) == 1)
    
        %         DS_sub6(k,n) = p_d*norm(reshape(D_sub6_sub6(k,k,n,:),[1,N_UE_sub6]))^2;
                DS_sub6(k,n) = p_d*(abs(D_sub6_sub6(k,k,n,n,n_idx)))^2;
                for nn = 1:N_UE_cell
        %             if (nn~=n)
        %             if (nn>n)
                    if (abs(D_sub6_sub6(k,k,nn,nn,n_idx))<abs(D_sub6_sub6(k,k,n,n,n_idx)))
        %                 MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,k,nn,:),[1,N_UE_sub6]))^2;
                        MSI_sub6(k,n) = MSI_sub6(k,n) + p_d*(abs(D_sub6_sub6(k,k,n,nn,n_idx)))^2;
                    end
                end
                for q = 1:K_mmW
                    if (q~=k && sub6ConnectionState(q)==1)
                      MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_mmW(k,q,n,:,n_idx),[1,N_UE_FWA]))^2;
                    end
                end
                for q = 1:K-K_mmW
                   MUI_sub6(k,n) = MUI_sub6(k,n) + p_d*norm(reshape(D_sub6_sub6(k,q,n,:,n_idx),[1,N_UE_cell]))^2;
                end
                snr_num_sub6(k,n) = DS_sub6(k,n);
                snr_den_sub6(k,n) = MSI_sub6(k,n) + MUI_sub6(k,n) + noise_sub6(k,n);
                % rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + scs(2)*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
                rate_dl(k+K_mmW) = rate_dl(k+K_mmW) + scs(n_idx)*TAU_FAC*log2(1+snr_num_sub6(k,n)/snr_den_sub6(k,n));
            end
        end
    end
end
end