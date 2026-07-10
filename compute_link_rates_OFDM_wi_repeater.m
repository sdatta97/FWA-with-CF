function rate_dl = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_dl_FWA, channel_interFWA)
M_sectors = params.numGNB;
K_FWA = params.numCPE;
K = M_sectors*params.numUE + params.numCPE;
BW = params.Band;
P_band = BW/M_sectors; %orthogonal band per sector (reuse factor = M_sectors)
TAU_FAC = params.preLogFactor;
N_BS = size(channel_dl,3);
N_UE = size(channel_dl,4);
N_CPE_FWA = size(channel_dl_FWA,4);
p_d = params.rho_tot; % 1*K;
D_FWA = params.D_FWA;
D_cell = params.D_cell;
BETA_interUE = params.BETA_interUE;
rep_gain = params.repeat_gain;
%Donor-side sector antenna of the network-controlled repeaters: amplitude
%gain of CPE q's best-pointing sector panel towards sector BS m
rep_donor_amp = sqrt(params.rep_donor_gain);
%BS-to-CPE large-scale gain (linear), from the shared gain matrix; used
%only for the repeater-selection metric
BETA_FWA = params.gainOverNoise_lin(:,1:K_FWA).*D_FWA;
if params.HW_IMPAIRMENTS
    Kt = params.Kt;
    Kr = params.Kr;
    K_rep_hw = params.Kt_rep*params.Kr_rep; %repeater AF-chain quality factor
else
    Kt = 1; Kr = 1; K_rep_hw = 1;
end
%Cellular UEs scheduled by each sector; the per-sector band is shared
%equally among them (OFDMA), so each UE gets P_band/sectorLoad
sectorLoad = sum(D_cell,2);
%Serving sector of each UE
Serv = cell(K-K_FWA,1);
for k = 1:K-K_FWA
    Serv{k} = find(D_cell(:,k)==1);
end
K_rep = params.num_repeater_per_cpe;
%Repeaters (CPEs) associated with each UE
Rep = cell(K-K_FWA,1);
for k = 1:K-K_FWA
    %repeater selection metric includes both repeater antenna gains:
    %donor (in rep_donor_gain) and service (in BETA_interUE)
    v = BETA_FWA(Serv{k},:)'.*params.rep_donor_gain(Serv{k},:)'.*BETA_interUE(1:K_FWA,k+K_FWA);
    %only CPEs enabled as repeaters (params.set_repeat) may be selected;
    %each UE associates with at most K_rep of them
    v(setdiff(1:K_FWA,params.set_repeat)) = -Inf;
    [vmax,servingCPEs] =  maxk(v,K_rep);
    Rep{k} = servingCPEs(isfinite(vmax));
end

%% Computing rates
DS_dl = zeros(K-K_FWA,N_UE);
MSI_dl = zeros(K-K_FWA,N_UE);
HI_dl = zeros(K-K_FWA,N_UE);
noise_dl = abs(sqrt(0.5)*(randn(K-K_FWA,N_UE) + 1j*randn(K-K_FWA,N_UE))).^2;
rate_dl = zeros(K-K_FWA,1);
for k = 1:K-K_FWA
    for n = 1:N_UE
        %receive channel of UE k split into direct and repeater-path
        %components: the useful signal keeps only sqrt(K_rep_hw) of the
        %repeater branch, the rest becomes repeater hardware distortion
        eff_dir = reshape(channel_dl(Serv{k},k,:,n),N_BS,1);
        eff_rep = zeros(N_BS,1);
        for kk = 1:numel(Rep{k})
            rep_idx = Rep{k}(kk);
            eff_rep = eff_rep + reshape(channel_dl_FWA(Serv{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(Serv{k},rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,n),[N_CPE_FWA,1]);
        end
        eff_channel = eff_dir + eff_rep;         %composite (beamforming basis)
        eff_sig = eff_dir + sqrt(K_rep_hw)*eff_rep; %coherent part after repeater HI
        ds_base = p_d*abs(eff_sig'*eff_channel./norm(eff_channel))^2;
        DS_dl(k,n) = Kr*Kt*ds_base;
        HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*ds_base ...
                   + (1-K_rep_hw)*p_d*abs(eff_rep'*eff_channel./norm(eff_channel))^2;
        for nn = 1:N_UE
            if (nn~=n)
                nn_eff_channel = reshape(channel_dl(Serv{k},k,:,nn),N_BS,1);
                 for kk = 1:numel(Rep{k})
                    rep_idx = Rep{k}(kk);
                    nn_eff_channel = nn_eff_channel + reshape(channel_dl_FWA(Serv{k},rep_idx,:,:),[N_BS,N_CPE_FWA])*rep_donor_amp(Serv{k},rep_idx)*rep_gain*reshape(channel_interFWA(rep_idx,k+K_FWA,:,nn),[N_CPE_FWA,1]);
                end
                if (norm(nn_eff_channel,'fro') < norm(eff_channel,'fro'))
                    msi_base = p_d*abs(eff_sig'*nn_eff_channel./norm(nn_eff_channel))^2;
                    MSI_dl(k,n) = Kr*Kt*msi_base;
                    HI_dl(k,n) = HI_dl(k,n) + (1-Kr*Kt)*msi_base ...
                               + (1-K_rep_hw)*p_d*abs(eff_rep'*nn_eff_channel./norm(nn_eff_channel))^2;
                end
            end
        end
        rate_dl(k) = rate_dl(k) + (P_band/sectorLoad(Serv{k}))*TAU_FAC*log2(1+DS_dl(k,n)/(MSI_dl(k,n)+HI_dl(k,n)+noise_dl(k,n)));
    end
end
end
