function rate_dl = compute_link_rates_OFDM(params,channel_dl)
M = params.numGNB;
K_FWA = params.numCPE;
K = M*params.numUE + params.numCPE;
BW = params.Band;
alpha = params.UE_split;
I_band = alpha*BW;
P_band = (1-alpha)*BW/M;
K_P = floor((1-alpha)*params.numUE);
K_I = K-K_FWA-M*K_P;
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D = params.D;
BETA = params.BETA;
BETA = BETA.*D;
if params.HW_IMPAIRMENTS
    Kt = params.Kt;
    Kr = params.Kr;
else
    Kt = 1; Kr = 1;
end
P_idxs = cell(M,1);
I_idxs = [];
for m = 1:M
    P_idxs{m,1} = find(BETA(m,:));
    I_idxs = [I_idxs,setdiff(find(D(m,:)),P_idxs{m,1})];
end
%Prepare cell to store the AP indices serving a specfic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specfic UE
NoServ = cell(K,1);
%Construc the above array and cells
for k = 1:K
    servingBSs = find(D(:,k)==1);
    NoservingBSs = find(D(:,k)==0);
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
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
HI_dl = zeros(K-K_FWA,N_UE);
noise_dl = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    for n = 1:N_UE
        if params.BEAM
            ds_base = p_d*abs(reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,n),N_BS,1)./norm(channel_dl(Serv{k},k,:,n),'fro'))^2;
            DS_dl(k,n) = Kr*Kt*ds_base;
            HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base;
            for nn = 1:N_UE
                if (nn~=n)
                    if (reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)<reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,n),N_BS,1))
                        msi_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)./norm(channel_dl(Serv{k},k,:,nn),'fro'))^2;
                        MSI_dl(k,n) = Kr*Kt*msi_base;
                        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base;
                    end
                    if ismember(k,I_idxs)
                        mci_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*reshape(channel_dl(NoServ{k},k,:,nn),N_BS,1)./norm(channel_dl(NoServ{k},k,:,nn),'fro'))^2;
                        MCI_dl(k,n) = Kr*Kt*mci_base;
                        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base;
                    end
                end
            end
            if ismember(k,I_idxs)
                mci_base2 = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*reshape(channel_dl(NoServ{k},k,:,n),N_BS,1)./norm(channel_dl(NoServ{k},k,:,n),'fro'))^2;
                MCI_dl(k,n) = MCI_dl(k,n) + Kr*Kt*mci_base2;
                HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base2;
            end
            if ismember(k,I_idxs)
                rate_dl(k) = rate_dl(k) + (I_band/(numel(I_idxs)/M))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
            else
<<<<<<< HEAD
                rate_dl(k) = rate_dl(k) + (P_band/(numel(P_idxs)/M))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
=======
                rate_dl(k) = rate_dl(k) + (P_band/numel(P_idxs(Serv{k},:)))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+noise_dl(k,n)));
>>>>>>> 0eee83c2230aaf6f881290ee81494ec12bb01b5a
            end
        else
            ds_base = p_d*abs(reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*ones(N_BS,1)./sqrt(N_BS))^2;
            DS_dl(k,n) = Kr*Kt*ds_base;
            HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base;
            for nn = 1:N_UE
                if (nn~=n)
                    if (reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)<reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,n),N_BS,1))
                        msi_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
                        MSI_dl(k,n) = Kr*Kt*msi_base;
                        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base;
                    end
                    if ismember(k,I_idxs)
                        mci_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
                        MCI_dl(k,n) = Kr*Kt*mci_base;
                        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base;
                    end
                end
            end
            if ismember(k,I_idxs)
<<<<<<< HEAD
                mci_base2 = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
                MCI_dl(k,n) = MCI_dl(k,n) + Kr*Kt*mci_base2;
                HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base2;
=======
                MCI_dl(k,n) = MCI_dl(k,n) + p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
>>>>>>> 0eee83c2230aaf6f881290ee81494ec12bb01b5a
            end
            if ismember(k,I_idxs)
                rate_dl(k) = rate_dl(k) + (I_band/(numel(I_idxs)/M))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
            else
<<<<<<< HEAD
                rate_dl(k) = rate_dl(k) + (P_band/(numel(P_idxs)/M))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
=======
                rate_dl(k) = rate_dl(k) + (P_band/numel(P_idxs(Serv{k},:)))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+noise_dl(k,n)));
>>>>>>> 0eee83c2230aaf6f881290ee81494ec12bb01b5a
            end
        end
    end
end
end