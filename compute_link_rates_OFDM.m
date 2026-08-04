function rate_dl = compute_link_rates_OFDM(params,channel_dl,channel_est_dl)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor; %scalar, or per-UE vector (per-class CSI overheads)
if isscalar(TAU_FAC)
    TAU_FAC = TAU_FAC*ones(K - params.numCPE,1);
end
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
%(matched filter per scheduled UE), built on
%the ESTIMATED channels (imperfect CSI when params.IMPERFECT_CSI is set);
%physical reception uses the true channels. The same matrices provide the
%victim's own stream beams and the beamformed inter-sector interference,
%averaged over the interferer's scheduled set (round-robin OFDMA)
W = cell(M_sectors,1);
for m = 1:M_sectors
    W{m} = beamMatrix(channel_est_dl,m,U{m},N_BS,N_UE);
end

%% Computing rates: MIMO log-det per UE
%Shannon capacity of UE k's N_UE x N_UE effective channel A (receive
%antennas x streams, streams beamformed per antenna) under colored
%interference-plus-noise:
%  rate_k = share*prelog*log2 det(I + Kt*Kr*A' R_k^{-1} A),
%R_k = beamformed inter-sector interference as a spatial covariance
%(averaged over the interferer's scheduled beams, round-robin OFDMA)
%+ diagonal EVM distortion + noise.
noise_dl = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    m = Serv{k};
    share = BW/sectorLoad(m);
    W_own = W{m}; U_own = U{m};
    int_sectors = setdiff(1:M_sectors,m); %reuse 1: every other sector interferes
    H_own = reshape(channel_dl(m,k,:,1:N_UE),[N_BS,N_UE]); %columns = receive antennas
    Wk = zeros(N_BS,N_UE); %the serving sector's beams for k's streams
    for n = 1:N_UE
        Wk(:,n) = getBeam(W_own,U_own,k,n,N_UE);
    end
    A = sqrt(p_d)*(H_own'*Wk); %A(nr,ns): stream ns as seen at antenna nr
    %inter-sector (co-channel) interference with the interferers'
    %actual beamforming towards their own scheduled UEs
    R_mci = zeros(N_UE);
    for mm = int_sectors
        if ~isempty(W{mm})
            H_int = reshape(channel_dl(mm,k,:,1:N_UE),[N_BS,N_UE]);
            B = sqrt(p_d)*(H_int'*W{mm});
            R_mci = R_mci + (B*B')/size(W{mm},2);
        end
    end
    R = diag(noise_dl(k,:)) + Kr*Kt*R_mci ...
      + diag((1-Kr*Kt)*(sum(abs(A).^2,2) + real(diag(R_mci))));
    rate_dl(k) = share*TAU_FAC(k)*logdet2(Kr*Kt*(A'*(R\A)));
end
end

function W = beamMatrix(channel_dl,m,U,N_BS,N_UE)
%Columns are the unit-norm matched-filter beams sector m applies towards
%its scheduled UEs (one per UE antenna stream)
if isempty(U)
    W = zeros(N_BS,0);
    return
end
W = zeros(N_BS,numel(U)*N_UE);
c = 0;
for q = U(:)'
    for nn = 1:N_UE
        c = c + 1;
        h = reshape(channel_dl(m,q,:,nn),N_BS,1);
        W(:,c) = h./norm(h);
    end
end
end

function w = getBeam(W,U,q,nn,N_UE)
%Beam a sector applies towards antenna nn of its scheduled UE q: a column
%of its beam matrix
if size(W,2) == 1
    w = W(:,1);
else
    w = W(:,(find(U==q,1)-1)*N_UE+nn);
end
end
