function rate_dl = compute_link_rates_MIMO_mmse(params,channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band; %frequency reuse 1: every sector uses the full band
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_CPE_FWA = size(channel_dl_FWA,4);
N_UE = size(channel_dl,4);
p_d = params.rho_tot; %total sector power. Constant-PSD convention: channels
%are normalized by FULL-band noise, so p_d*|h|^2 is the PSD ratio - exact
%on any subband; the bandwidth share enters only through the rate prelog
D_FWA = params.D_FWA;
set_repeat = params.set_repeat;
SI_cancel_factor = params.SI_cancel_factor;
if params.HW_IMPAIRMENTS
    Kt     = params.Kt;
    Kr_FWA  = params.Kr;
    Kr_cell = params.Kr;
else
    Kt = 1; Kr_FWA = 1; Kr_cell = 1;
end

%Serving sector(s) of each user
Serv = cell(K,1);
for k = 1:K
    Serv{k} = find(D_FWA(:,k)==1);
end

%% initialization of c
D_FWA_FWA = zeros(K_FWA,K_FWA,N_CPE_FWA,N_CPE_FWA);
D_FWA_Cell = zeros(K_FWA,K-K_FWA,N_CPE_FWA,N_UE);
D_Cell_FWA = zeros(K-K_FWA,K_FWA,N_UE,N_CPE_FWA);
D_Cell_Cell = zeros(K-K_FWA,K-K_FWA,N_UE,N_UE);
dl_mmse_precoder_FWA = zeros(size(channel_est_dl_FWA));
dl_mmse_precoder = zeros(size(channel_est_dl));
scaling_LP_mmse = zeros(M_sectors,K);
%MMSE precoding, structured like the OFDM beam matrices: the MMSE filter
%basis inv_matrix is identical for every user a sector precodes, so it
%is built once per sector and applied to each served user's channel.
%Beamforming operates on the ESTIMATED channels (imperfect CSI when
%params.IMPERFECT_CSI is set); physical reception below uses the true ones
for m = 1:M_sectors
    inv_matrix = eye(N_BS);
    for q = 1:K_FWA
        if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])*reshape(channel_est_dl_FWA(m,q,:,:),[N_BS,N_CPE_FWA])';
        end
    end
    for q = 1:K-K_FWA
        if ismember(m,Serv{q+K_FWA})
            inv_matrix = inv_matrix + p_d*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])*reshape(channel_est_dl(m,q,:,:),[N_BS,N_UE])';
        end
    end
    for k = 1:K_FWA
        if ~ismember(k,set_repeat)
            dl_mmse_precoder_FWA(m,k,:,:) = p_d*inv_matrix\(reshape(channel_est_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA]));
            if ismember(m,Serv{k})
                scaling_LP_mmse(m,k) = norm(dl_mmse_precoder_FWA(m,k,:,:),'fro')^2;
            end
        end
    end
    for k = 1:K-K_FWA
        dl_mmse_precoder(m,k,:,:) = p_d*inv_matrix\(reshape(channel_est_dl(m,k,:,:),[N_BS,N_UE]));
        if ismember(m,Serv{k+K_FWA})
            scaling_LP_mmse(m,k+K_FWA) = norm(dl_mmse_precoder(m,k,:,:),'fro')^2;
        end
    end
end
for m = 1:M_sectors
    for k = 1:K_FWA
        if ismember(m,Serv{k}) && ~ismember(k,set_repeat)
            dl_mmse_precoder_FWA(m,k,:,:) = reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])./sqrt(scaling_LP_mmse(m,k));
        end
    end
    for k = 1:K-K_FWA
        if ismember(m,Serv{k+K_FWA})
            dl_mmse_precoder(m,k,:,:) = reshape(dl_mmse_precoder(m,k,:,:),[N_BS,N_UE])./sqrt(scaling_LP_mmse(m,k+K_FWA));
        end
    end
end
eta_eq = zeros(M_sectors,K);
for m = 1:M_sectors
    term = 0;
    for k = 1:K
        if ismember(m,Serv{k})
            if (k<=K_FWA)
                if ~ismember(k,set_repeat)
                    term = term + trace(reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])*reshape(dl_mmse_precoder_FWA(m,k,:,:),[N_BS,N_CPE_FWA])');
                end
            else
                term = term + trace(reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])*reshape(dl_mmse_precoder(m,k-K_FWA,:,:),[N_BS,N_UE])');
            end
        end
    end
    if (term > 0)
        setA = 1:K;
        eta_eq(m,setA(~ismember(setA,set_repeat))) = (1/term)*D_FWA(m,setA(~ismember(setA,set_repeat)));
    end
end
%Effective gains: the channels channel_dl_FWA/channel_dl already carry
%the sector antenna pattern of BS m towards each user (via the
%pattern-inclusive large-scale gains inside R_gNB), so both the desired
%links and every inter-sector interference term below are sectored
for k = 1:K_FWA
    if ~ismember(k,set_repeat)
        for q = 1:K_FWA
            for m = 1:M_sectors
                if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
                    D_FWA_FWA(k,q,:,:) = reshape(D_FWA_FWA(k,q,:,:),[N_CPE_FWA,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
                end
            end
        end
        for q = 1:K-K_FWA
            for m = 1:M_sectors
                if ismember(m,Serv{q+K_FWA})
                    D_FWA_Cell(k,q,:,:) = reshape(D_FWA_Cell(k,q,:,:),[N_CPE_FWA,N_UE]) + sqrt(p_d*eta_eq(m,q+K_FWA))*reshape(channel_dl_FWA(m,k,:,:),[N_BS,N_CPE_FWA])'*reshape(dl_mmse_precoder(m,q,:,:),[N_BS,N_UE]);
                end
            end
        end
    end
end
for k = 1:K-K_FWA
    for q = 1:K_FWA
        for m = 1:M_sectors
            if ismember(m,Serv{q}) && ~ismember(q,set_repeat)
                D_Cell_FWA(k,q,:,:) = reshape(D_Cell_FWA(k,q,:,:),[N_UE,N_CPE_FWA]) + sqrt(p_d*eta_eq(m,q))*reshape(channel_dl(m,k,:,:),[N_BS,N_UE])'*reshape(dl_mmse_precoder_FWA(m,q,:,:),[N_BS,N_CPE_FWA]);
            end
        end
    end
    for q = 1:K-K_FWA
        for m = 1:M_sectors
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
    if ~ismember(k,set_repeat)
        for n = 1:N_CPE_FWA
            DS_FWA(k,n)  = Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,n)))^2;
            HI_FWA(k,n)  = (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,n)))^2;
            for nn = 1:N_CPE_FWA
                if (abs(D_FWA_FWA(k,k,nn,nn))<abs(D_FWA_FWA(k,k,n,n)))
                    MSI_FWA(k,n) = MSI_FWA(k,n) + Kr_FWA*Kt*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                    HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*(abs(D_FWA_FWA(k,k,n,nn)))^2;
                end
            end
            for q = 1:K_FWA
                if (q~=k) && ~ismember(q,set_repeat)
                    MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*SI_cancel_factor*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                    HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*SI_cancel_factor*norm(reshape(D_FWA_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                end
            end
            for q = 1:K-K_FWA
                MUI_FWA(k,n) = MUI_FWA(k,n) + Kr_FWA*Kt*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
                HI_FWA(k,n)  = HI_FWA(k,n)  + (1-Kr_FWA*Kt)*norm(reshape(D_FWA_Cell(k,q,n,:),[1,N_UE]))^2;
            end
            snr_num_FWA(k,n) = DS_FWA(k,n);
            snr_den_FWA(k,n) = MSI_FWA(k,n) + MUI_FWA(k,n) + HI_FWA(k,n) + noise_FWA(k,n);
            rate_dl(k) = rate_dl(k) + BW*TAU_FAC*log2(1+snr_num_FWA(k,n)/snr_den_FWA(k,n));
        end
    end
end
for k = 1:K-K_FWA
    for n = 1:N_UE
        DS_Cell(k,n)  = Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,n)))^2;
        HI_Cell(k,n)  = (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,n)))^2;
        for nn = 1:N_UE
            if (abs(D_Cell_Cell(k,k,nn,nn))<abs(D_Cell_Cell(k,k,n,n)))
                MSI_Cell(k,n) = MSI_Cell(k,n) + Kr_cell*Kt*(abs(D_Cell_Cell(k,k,n,nn)))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*(abs(D_Cell_Cell(k,k,n,nn)))^2;
            end
        end
        for q = 1:K_FWA
            if ~ismember(q,set_repeat)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*norm(reshape(D_Cell_FWA(k,q,n,:),[1,N_CPE_FWA]))^2;
            end
        end
        for q = 1:K-K_FWA
            if (q~=k)
                MUI_Cell(k,n) = MUI_Cell(k,n) + Kr_cell*Kt*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
                HI_Cell(k,n)  = HI_Cell(k,n)  + (1-Kr_cell*Kt)*norm(reshape(D_Cell_Cell(k,q,n,:),[1,N_UE]))^2;
            end
        end
        snr_num_Cell(k,n) = DS_Cell(k,n);
        snr_den_Cell(k,n) = MSI_Cell(k,n) + MUI_Cell(k,n) + HI_Cell(k,n) + noise_Cell(k,n);
        rate_dl(k+K_FWA) = rate_dl(k+K_FWA) + BW*TAU_FAC*log2(1+snr_num_Cell(k,n)/snr_den_Cell(k,n));
    end
end
end
