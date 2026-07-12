function rate_dl = compute_link_rates_OFDM(params,channel_dl,channel_est_dl)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; %total sector power. Constant-PSD convention: channels
%are normalized by FULL-band noise, so p_d*|h|^2 is the PSD ratio - exact
%on any subband; the bandwidth share enters only through the rate prelog
D = params.D;
if params.HW_IMPAIRMENTS
    Kt = params.Kt;
    Kr = params.Kr;
else
    Kt = 1; Kr = 1;
end
%Serving sector of each UE and per-sector scheduled sets; each sector
%shares the full band equally among its UEs (OFDMA)
Serv = cell(K,1);
for k = 1:K
    Serv{k} = find(D(:,k)==1);
end
U = cell(M_sectors,1);
for m = 1:M_sectors
    U{m} = find(D(m,:)==1);
end
sectorLoad = sum(D,2);
%Beamforming vectors each sector applies towards its own scheduled UEs
%(matched filter for BEAM=1, fixed broadside beam for BEAM=0), built on
%the ESTIMATED channels (imperfect CSI when params.IMPERFECT_CSI is set);
%physical reception uses the true channels. The same matrices provide the
%victim's own DS/MSI beams and the beamformed inter-sector interference,
%averaged over the interferer's scheduled set (round-robin OFDMA)
W = cell(M_sectors,1);
for m = 1:M_sectors
    W{m} = beamMatrix(channel_est_dl,m,U{m},N_BS,N_UE,params.BEAM);
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
        h_n = reshape(channel_dl(m,k,:,n),N_BS,1);
        ds_base = p_d*abs(h_n'*getBeam(W_own,U_own,k,n,N_UE))^2;
        DS_dl(k,n) = Kr*Kt*ds_base;
        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base;
        for nn = 1:N_UE
            if (nn~=n)
                h_nn = reshape(channel_dl(m,k,:,nn),N_BS,1);
                if (h_nn'*h_nn < h_n'*h_n)
                    msi_base = p_d*abs(h_n'*getBeam(W_own,U_own,k,nn,N_UE))^2;
                    MSI_dl(k,n) = Kr*Kt*msi_base;
                    HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base;
                end
            end
        end
        %inter-sector (co-channel) interference with the interferers'
        %actual beamforming towards their own scheduled UEs
        mci_base = 0;
        for mm = int_sectors
            if ~isempty(W{mm})
                h_int = reshape(channel_dl(mm,k,:,n),N_BS,1);
                mci_base = mci_base + p_d*mean(abs(h_int'*W{mm}).^2);
            end
        end
        MCI_dl(k,n) = Kr*Kt*mci_base;
        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*mci_base;
        rate_dl(k) = rate_dl(k) + share*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+MCI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
    end
end
end

function W = beamMatrix(channel_dl,m,U,N_BS,N_UE,BEAM)
%Columns are the unit-norm beams sector m applies towards its scheduled
%UEs (one per UE antenna stream)
if isempty(U)
    W = zeros(N_BS,0);
elseif BEAM
    W = zeros(N_BS,numel(U)*N_UE);
    c = 0;
    for q = U(:)'
        for nn = 1:N_UE
            h = reshape(channel_dl(m,q,:,nn),N_BS,1);
            c = c + 1;
            W(:,c) = h./norm(h);
        end
    end
else
    W = ones(N_BS,1)./sqrt(N_BS); %fixed beam, identical for every UE
end
end

function w = getBeam(W,U,q,nn,N_UE)
%Beam a sector applies towards antenna nn of its scheduled UE q: a column
%of its beam matrix (the single shared column under the fixed-beam model)
if size(W,2) == 1
    w = W(:,1);
else
    w = W(:,(find(U==q,1)-1)*N_UE+nn);
end
end
