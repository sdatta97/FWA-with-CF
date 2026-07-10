function rate_dl = compute_link_rates_OFDM(params,channel_dl)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band;
P_band = BW/M_sectors; %orthogonal band per sector (reuse factor = M_sectors)
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; % 1*K;
D = params.D;
if params.HW_IMPAIRMENTS
    Kt = params.Kt;
    Kr = params.Kr;
else
    Kt = 1; Kr = 1;
end
%Cellular UEs scheduled by each sector; the per-sector band is shared
%equally among them (OFDMA), so each UE gets P_band/sectorLoad
sectorLoad = sum(D,2);
%Serving sector of each UE
Serv = cell(K,1);
for k = 1:K
    Serv{k} = find(D(:,k)==1);
end

%% Computing rates
DS_dl = zeros(K-K_FWA,N_UE);
MSI_dl = zeros(K-K_FWA,N_UE);
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
                if (nn~=n) && (reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)<reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,n),N_BS,1))
                    msi_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)./norm(channel_dl(Serv{k},k,:,nn),'fro'))^2;
                    MSI_dl(k,n) = Kr*Kt*msi_base;
                    HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base;
                end
            end
        else
            ds_base = p_d*abs(reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*ones(N_BS,1)./sqrt(N_BS))^2;
            DS_dl(k,n) = Kr*Kt*ds_base;
            HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base;
            for nn = 1:N_UE
                if (nn~=n) && (reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,nn),N_BS,1)<reshape(channel_dl(Serv{k},k,:,n),N_BS,1)'*reshape(channel_dl(Serv{k},k,:,n),N_BS,1))
                    msi_base = p_d*abs((reshape(channel_dl(Serv{k},k,:,n),N_BS,1))'*ones(N_BS,1)./sqrt(N_BS))^2;
                    MSI_dl(k,n) = Kr*Kt*msi_base;
                    HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base;
                end
            end
        end
        rate_dl(k) = rate_dl(k) + (P_band/sectorLoad(Serv{k}))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
    end
end
end
