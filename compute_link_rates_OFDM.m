function rate_dl = compute_link_rates_OFDM(params,channel_dl)
M = params.numGNB;
K_FWA = params.numCPE;
K = M*params.numUE + params.numCPE;
BW = params.Band;
alpha = params.UE_split;
I_band = alpha*BW;
P_band = (1-alpha)*BW/M;
K_P = floor((1-alpha)*params.numUE);
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
    [~,P_idxs{m,1}] = maxk(BETA(m,:),K_P);
    I_idxs = [I_idxs,setdiff(find(D(m,:)),P_idxs{m,1})];
end
%Prepare cell to store the AP indices serving a specfic UE
Serv = cell(K,1);
%Prepare cell to store the AP indices not serving a specfic UE
NoServ = cell(K,1);
%Construct the above array and cells
for k = 1:K
    servingBSs = find(D(:,k)==1);
    NoservingBSs = find(D(:,k)==0);
    
    Serv{k} = servingBSs;
    NoServ{k} = NoservingBSs;
end

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
                rate_dl(k) = rate_dl(k) + (P_band/numel(P_idxs(Serv{k},:)))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
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
                mci_base2 = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
                MCI_dl(k,n) = MCI_dl(k,n) + Kr*Kt*mci_base2;
                HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base2;
            end
            if ismember(k,I_idxs)
                rate_dl(k) = rate_dl(k) + (I_band/(numel(I_idxs)/M))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
            else
                rate_dl(k) = rate_dl(k) + (P_band/numel(P_idxs(Serv{k},:)))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
            end
        end
    end
end
end