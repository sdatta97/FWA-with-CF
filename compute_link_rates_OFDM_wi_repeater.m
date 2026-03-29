function rate_dl = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_dl_FWA, channel_interFWA)
M = params.numGNB;
K_FWA = params.numCPE;
K = M*params.numUE + params.numCPE;
BW = params.Band;
alpha = 0;
I_band = alpha*BW;
P_band = (1-alpha)*BW/M;
K_P = floor((1-alpha)*params.numUE);
K_I = K-K_FWA-M*K_P;
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
N_CPE_FWA = size(channel_dl_FWA,4);
p_d = params.rho_tot; % 1*K;
D = params.D;
BETA = params.BETA;
BETA_interUE = params.BETA_interUE;
rep_gain = params.repeat_gain;
BETA = BETA.*D;
P_idxs = cell(M,1);
I_idxs = [];
for m = 1:M
    P_idxs{m,1} = find(BETA(m,1+K_FWA:end));
    % [~,P_idxs_tmp] = mink(BETA(m,K_FWA+non_zero_idxs),min(K_P/M,numel(non_zero_idxs)));
    I_idxs = [I_idxs,setdiff(find(D(m,1+K_FWA:end)),P_idxs{m,1})];
end
%Prepare cell to store the AP indices serving a specfic UE
Serv = cell(K-K_FWA,1);
%Prepare cell to store the AP indices not serving a specfic UE
NoServ = cell(K-K_FWA,1);
%Construc the above array and cells
for k = 1:K-K_FWA
    servingBSs = find(D(:,k+K_FWA)==1);
    NoservingBSs = find(D(:,k+K_FWA)==0);
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
end
K_rep = params.num_repeater_per_cpe;
%Prepare cell to store the CPE indices serving a specfic UE
Rep = cell(K-K_FWA,1);
%Prepare cell to store the CPE indices not serving a specfic UE
NoRep = cell(K-K_FWA,1);
%Construc the above array and cells
for k = 1:K-K_FWA
    % v1 = BETA(Serv{k},1:K_FWA)'.*BETA_interUE(1:K_FWA,k+K_FWA);
    % v2 = BETA(Serv{k},1:K_FWA)'.*BETA_interUE(1:K_FWA,k+K_FWA) - BETA(NoServ{k},1:K_FWA)'.*BETA_interUE(1:K_FWA,k+K_FWA);
    % [~,servingCPEs1] =  maxk(v1,K_rep);
    % [~,servingCPEs2] =  maxk(v2,K_rep);
    % if (servingCPEs2 ~= servingCPEs1)
    %     disp("hi");
    % else
    %     servingCPEs = servingCPEs1;
    % end
    v = BETA(Serv{k},1:K_FWA)'.*BETA_interUE(1:K_FWA,k+K_FWA);
    [~,servingCPEs] =  maxk(v,K_rep);
    NoservingCPEs = setdiff(params.set_repeat,servingCPEs);
    Rep{k} = servingCPEs;
    NoRep{k} = NoservingCPEs; 
end

% %% initialization of c
% D_Cell_Cell = zeros(K-K_FWA,K-K_FWA,N_UE,N_UE);
% dl_mmse_precoder = zeros(size(channel_est_dl));
% scaling_LP_mmse = zeros(M,K-K_FWA);
% for m = 1:M
%     for k = 1:K-K_FWA
%         % inv_matrix = noiseVariance*eye(Ntx);
%         inv_matrix = eye(N_BS);
%         for q = 1:K-K_FWA
%             if ismember(m,Serv{q+K_FWA})
%                 inv_matrix = inv_matrix +  p_d*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_dl(m,q,:,:),[N_BS,N_UE])';
%             end
%         end
%         dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE]) + p_d*inv_matrix\(reshape(channel_dl(m,k,:,:),[N_BS,N_UE]));
%         if ismember(m,Serv{k+K_FWA})
%             scaling_LP_mmse(m,k+K_FWA) = scaling_LP_mmse(m,k+K_FWA) + norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
%         end
%     end
% end
% for m = 1:M
%     for k = 1:K-K_FWA
%         if ismember(m,Serv{k+K_FWA})
%             dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE])./sqrt(scaling_LP_mmse(m,k+K_FWA));
%         end
%     end
% end
% eta_eq = zeros(M,K-K_FWA);
% for m = 1:M
%     term = 0;
%     for k = 1:K-K_FWA
%         if ismember(m,Serv{k+K_FWA})
%             term = term + trace(reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])*reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])');
%         end
%     end
%     if (term > 0)
%         eta_eq(m,:) = (1/term)*D(m,:);
%     end
% end
% for k = 1:K-K_FWA
%     for q = 1:K-K_FWA
%         for m = 1:M
%             if ismember(m,Serv{q+K_FWA})
%                 D_Cell_Cell(k,q,:,:) = reshape(D_Cell_Cell(k,q,:,:),[N_UE,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
%             end
%         end
%     end
% end
%% Computing rates
DS_dl = zeros(K-K_FWA,N_UE);
MSI_dl = zeros(K-K_FWA,N_UE);
MCI_dl = zeros(K-K_FWA,N_UE); 
noise_dl = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    for n = 1:N_UE
        eff_channel = reshape(channel_dl(Serv{k},k,:,n),N_BS,1);
        for kk = 1:numel(Rep{k})
            rep_idx = Rep{k}(kk);
            eff_channel = eff_channel + reshape(channel_dl_FWA(Serv{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[N_CPE_FWA,1]);
        end
        DS_dl(k,n) = p_d*abs(eff_channel'*eff_channel./norm(eff_channel))^2;
        for nn = 1:N_UE
            if (nn~=n)
                nn_eff_channel = reshape(channel_dl(Serv{k},k,:,nn),N_BS,1);
                 for kk = 1:numel(Rep{k})
                    rep_idx = Rep{k}(kk);
                    nn_eff_channel = nn_eff_channel + reshape(channel_dl_FWA(Serv{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,nn),[N_CPE_FWA,1]);
                end
                if (norm(nn_eff_channel,'fro') < norm(eff_channel,'fro'))
                    MSI_dl(k,n) = p_d*abs(eff_channel'*nn_eff_channel./norm(nn_eff_channel))^2;
                end
                if ismember(k,I_idxs)
                    mci_eff_channel = reshape(channel_dl(NoServ{k},k,:,nn),N_BS,1);
                    for kk = 1:numel(Rep{k})
                        rep_idx = Rep{k}(kk);
                        mci_eff_channel = mci_eff_channel + reshape(channel_dl_FWA(NoServ{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,nn),[N_CPE_FWA,1]);
                    end
                    MCI_dl(k,n) = p_d*abs(eff_channel'*mci_eff_channel./norm(mci_eff_channel))^2;
                end
            end
        end
        if ismember(k,I_idxs)   
            mci_eff_channel = reshape(channel_dl(NoServ{k},k,:,n),N_BS,1);
            for kk = 1:numel(Rep{k})
                rep_idx = Rep{k}(kk);
                mci_eff_channel = mci_eff_channel + reshape(channel_dl_FWA(NoServ{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[N_CPE_FWA,1]);
            end
            MCI_dl(k,n) = MCI_dl(k,n) + p_d*abs(eff_channel'*mci_eff_channel./norm(mci_eff_channel))^2;
        end
        if ismember(k,I_idxs)
            rate_dl(k) = rate_dl(k) + I_band*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+noise_dl(k,n)));
        else
            rate_dl(k) = rate_dl(k) + (P_band/numel(P_idxs(Serv{k},:)))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+noise_dl(k,n)));
        end
    end
end
end