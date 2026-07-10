close all;
clear;
tStart = tic;
%PackingAnalysis: per-user downlink rate DISTRIBUTIONS for cellular UEs
%vs FWA CPEs, quantifying how much closer to capacity each service can be
%safely packed. Capacity planning dimensions for the busy-hour worst case:
%operators trigger expansion at ~70% PRB utilization (Survey on 5G Network
%Expansion Methods, IUCC 2021,
%https://hpcn.exeter.ac.uk/iucc2021/proceedings/pdfs/IUCC-CIT-DSCI-SmartCNS2021-40WP54zLa9Wagib9WOs48p/666700a524/666700a524.pdf)
%and keep busy-hour PRB utilization below ~80-85%
%(https://www.telecomhall.net/t/capacity-guidelines-number-of-users-prb-utilization-data-volume-throughput/13458;
%https://www.rfwireless-world.com/terminology/4g-lte-prb-utilization).
%The margin exists because cellular demand is a RANDOM population of
%mobile users; FWA demand is measured per CPE at installation, and the
%carrier may decline service on poor links (FCC 24-27,
%https://docs.fcc.gov/public/attachments/FCC-24-27A1.pdf). The safe load
%fraction of a service is estimated as f(eps) = q_eps(rate)/mean(rate),
%which is invariant to the bandwidth split, so all rates here are
%computed on the full band. Cellular pools rates across seeds AND
%realizations (population uncertainty: aggregate with the other SLURM
%array outputs in dataProcess/plot_packing.m); FWA uses the per-drop
%temporal distribution only (installation is known). The ratio
%f_FWA/f_cell is the packing multiple of revenue potential.
%
%One SLURM array task = one deployment drop (seed), matching the
%SimulationMain.m parallelization structure.
aID = getenv('SLURM_ARRAY_TASK_ID');
if (isempty(aID))
    warning('aID is empty. Trying SLURM ID.')
    aID = getenv('SLURM_ARRAY_TASK_ID');
end
if(isempty(aID))
    warning('aID is empty. Replacing it with 0084.')
    aID = '0084';
end
%RNG seed.
rng(str2double(aID),'twister');
%% Mobility params (TR 38.901 Table 7.2-5)
params.MOBILE = 1;
params.ue_velocity_outdoor = 40/3.6; %outdoor (in-car) UTs: 40 km/h
params.ue_velocity_indoor = 3/3.6;   %indoor UTs: 3 km/h
params.Ts = 1.67e-5; %channel sample time (s)
params.fc = 3.7e9; %carrier frequency (Hz)
params.c = 3e8; %speed of light (m/s)
params.BEAM = 0;

numCPE_all = 50; %CPEs per site
Band = 100e6;
params.noiseFigure = 9; %receiver noise figure (dB)

%% Suburban Macro (SMa) deployment, 3GPP TR 38.901 Table 7.2-5
params.ISD = 1299;
params.deployRange = params.ISD/2;
params.complexRadius = 200;
params.coverageRange_sub6 = 2*params.deployRange + 30;
params.min_dist_2D = 35;
params.num_antennas_per_gNB = 64;
params.rho_tot = 10^(0.1*49); %49 dBm per sector (TR 38.901 clause 7.8)
params.N_UE_FWA = 8;
params.N_UE_cell = 2;
params.hr = 1.5;
params.hr_cpe = 4.5;
params.ht_bs = 35;
%% Sectored BS antennas (TR 38.901 Table 7.3-1)
params.sectors_per_site = 3;
sector_offsets = [30; 150; 270];
params.theta_3dB = 65;
params.phi_3dB = 65;
params.SLA_V = 30;
params.A_max = 30;
params.G_ant_max = 8;
if params.ISD <= 1500
    params.tilt_zenith = 95;
else
    params.tilt_zenith = 92;
end
%% SMa propagation (TR 38.901 Clause 7.4)
params.h_bldg = 10;
params.W_street = 10;
params.r_vegetation = 0.10;
params.indoor_UT_ratio = 0.8;
%Active-user densities, as in SimulationMain.m
lambda_UE_res = 150; %active UEs per km^2, residential block
lambda_UE_road = 4;  %active in-car UEs per km^2, roads outside
params.Lmax = 1;
params.preLogFactor = 1;
params.HW_IMPAIRMENTS = 1;
%EVM-based impairment factors (TS 38.104 Table 6.5.2-1; TS 38.101-1;
%repeater conformance TS 38.115-2), as in SimulationMain.m
params.Kt = 1 - 0.035^2;
params.Kr = 1 - 0.08^2;
params.Kt_rep = 1 - 0.08^2;
params.Kr_rep = 1 - 0.08^2;
params.IMPERFECT_CSI = 1;
params.csi_rs_offset_dB = 0;
%FWA residual interference after nulling, gamma_I (see SimulationMain.m
%for the Tarana / measured-null-depth citations)
SI_cancel_dB = -20;
%NCR amplification (TR 38.867: 90 dB gain, 33 dBm max output power)
params.repeat_gain = 90;
params.rep_max_pow_dBm = 33;
params.num_repeater_tot = 2; %enabled NCR pool for the cellular phase
params.num_repeater_per_cpe = 1;
params.rep_assist_frac = 0.2;
%Snapshots per drop: each realization is a time instant dt_snap apart
%along the cellular UEs' motion trajectories (see the mobility model
%below); FWA CPEs are stationary throughout
nbrOfRealizations = 50;
%% Temporal mobility model for the cellular UEs
%Straight-line constant-speed motion with a random heading per user (the
%Random Direction model: Camp, Boleng & Davies, "A survey of mobility
%models for ad hoc network research", Wireless Comm. and Mobile
%Computing, 2002, https://onlinelibrary.wiley.com/doi/10.1002/wcm.72),
%which is also the 3GPP / ITU-R M.2412 evaluation convention of straight
%trajectories at constant speed. Road UEs move at 40 km/h in the outer
%ring - the straight-segment core of the Freeway highway mobility model
%(Bai, Sadagopan & Helmy, "IMPORTANT: a framework to systematically
%analyze the impact of mobility on performance of routing protocols for
%adhoc networks", IEEE INFOCOM 2003,
%https://ieeexplore.ieee.org/document/1208920) - and indoor UEs walk at
%3 km/h inside the housing complex (speeds per TR 38.901 Table 7.2-5).
%Users reflect specularly at their region boundaries, so road UEs stay
%on the roads and indoor UEs stay in the complex. Across snapshots the
%shadow fading, LOS states and O2I losses evolve consistently with the
%traveled distance (TR 38.901 Clause 7.6.3 spatial consistency; see
%generateSetup.m), sector association follows the motion (ideal
%handover), and small-scale fading is independent (dt_snap far exceeds
%the channel coherence time). Within each snapshot the CSI-lag channel
%aging keeps the Jakes autocorrelation model.
params.dt_snap = 1; %time between snapshots (s)
params.TEMPORAL_MOBILITY = 1; %evolve UE positions/large-scale state across
                              %snapshots (gates generateSetup.m's evolveMode)
%Planning outage targets for the in-run summary
eps_grid = [0.01 0.02 0.05 0.1 0.2];

%% Deployment: same construction as SimulationMain.m
params.M_sites = 2;
M_sites = params.M_sites;
S = params.sectors_per_site;
params.numGNB = M_sites*S;
M_sectors = params.numGNB;
Rsite = [0; params.ISD];
angleSite = 2*pi*rand(M_sites,1);
siteLocations = [Rsite.*cos(angleSite), Rsite.*sin(angleSite)];
params.siteLocations = siteLocations;
params.locationsBS = kron(siteLocations, ones(S,1));
params.sector_boresights = repmat(sector_offsets, M_sites, 1);
numCPE_tot = M_sites*numCPE_all;
RCPE = sqrt(params.min_dist_2D^2 + (params.complexRadius^2 - params.min_dist_2D^2)*rand(numCPE_tot,1));
angleCPE = 2*pi*rand(numCPE_tot,1);
CPE_locations = kron(siteLocations, ones(numCPE_all,1)) + [RCPE.*cos(angleCPE), RCPE.*sin(angleCPE)];
A_complex = pi*((params.complexRadius/1000)^2 - (params.min_dist_2D/1000)^2);
A_road = pi*((params.deployRange/1000)^2 - (params.complexRadius/1000)^2);
numUE_in = ceil(lambda_UE_res*A_complex/S);
numUE_out = ceil(lambda_UE_road*A_road/S);
params.numUE = numUE_in + numUE_out;
numUE = params.numUE;
UE_locations = zeros(M_sectors*numUE,2);
isIndoorUE = false(M_sectors*numUE,1);
for idxSite = 1:M_sites
    Nin = S*numUE_in;
    Nout = S*numUE_out;
    Rin = sqrt(params.min_dist_2D^2 + (params.complexRadius^2 - params.min_dist_2D^2)*rand(Nin,1));
    Rout = sqrt(params.complexRadius^2 + (params.deployRange^2 - params.complexRadius^2)*rand(Nout,1));
    angleUE = 2*pi*rand(Nin+Nout,1);
    RUE = [Rin; Rout];
    base = (idxSite-1)*(Nin+Nout);
    UE_locations(base+(1:Nin+Nout),:) = siteLocations(idxSite,:) + [RUE.*cos(angleUE), RUE.*sin(angleUE)];
    isIndoorUE(base+(1:Nin)) = true;
end
params.UE_locations = UE_locations;
params.isIndoorUE = isIndoorUE;
params.numCPE = numCPE_tot;
params.CPE_locations = CPE_locations;
params.Band = Band;
[gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,~,D_FWA,D_cell,~,~,~,~,isIndoor,~,~,repDonorGaindB,mobState] = generateSetup(params,str2double(aID));
params.rep_donor_gain = db2pow(repDonorGaindB);
params.gainOverNoise_lin = db2pow(gainOverNoisedB);
noise_mW = 10^(0.1*(-174 + 10*log10(Band) + params.noiseFigure));
P_in_rep_mW = params.rho_tot*noise_mW*sum(params.gainOverNoise_lin(:,1:numCPE_tot).*params.rep_donor_gain,1)';
params.repeat_gain_amp = sqrt(min(10^(0.1*params.repeat_gain), ...
    10^(0.1*params.rep_max_pow_dBm)./P_in_rep_mW));
v_ue = params.ue_velocity_outdoor*ones(M_sectors*numUE,1);
v_ue(isIndoor(numCPE_tot+1:end)) = params.ue_velocity_indoor;
params.mob_rho = besselj(0,2*pi*params.Ts*v_ue*params.fc/params.c);
params.BETA_interUE = db2pow(gainOverNoisedB_ue);
%Enabled NCR pool: coverage-vote selection, as in SimulationMain.m
BETA_FWA_assoc = db2pow(gainOverNoisedB(:,1:numCPE_tot)).*D_FWA;
votes = zeros(numCPE_tot,1);
for kUE = 1:M_sectors*numUE
    m_serv = find(D_cell(:,kUE)==1,1);
    v_rep = BETA_FWA_assoc(m_serv,:)'.*params.rep_donor_gain(m_serv,:)'.*params.BETA_interUE(1:numCPE_tot,numCPE_tot+kUE);
    [v_best,best_rep] = max(v_rep);
    if isfinite(v_best) && v_best > 0
        votes(best_rep) = votes(best_rep) + 1;
    end
end
[~,vote_order] = sort(votes,'descend');
%NCR enablement is a semi-static network configuration: the pool is
%chosen once at the first snapshot and held over the motion timeline
set_repeat_pool = vote_order(1:min(params.num_repeater_tot,numCPE_tot));

%% Per-snapshot rates: at each time instant the cellular UEs have moved,
%the large-scale setup is re-evaluated consistently, and one small-scale
%realization is drawn per snapshot for all three populations. The safe
%load fraction f = q_eps/mean is bandwidth-scale-invariant, so all rates
%use the full band
K_FWA = numCPE_tot;
rate_cell_rep = zeros(M_sectors*numUE,nbrOfRealizations);
rate_cell_plain = zeros(M_sectors*numUE,nbrOfRealizations);
rate_FWA = zeros(K_FWA,nbrOfRealizations);
UEpos_t = UE_locations(:,1) + 1i*UE_locations(:,2);
siteCenterUE = repelem(siteLocations(:,1) + 1i*siteLocations(:,2), S*numUE, 1);
v_move = params.ue_velocity_outdoor*ones(M_sectors*numUE,1);
v_move(isIndoorUE) = params.ue_velocity_indoor;
heading = 2*pi*rand(M_sectors*numUE,1);
r_lo = params.min_dist_2D*ones(M_sectors*numUE,1); r_lo(~isIndoorUE) = params.complexRadius;
r_hi = params.complexRadius*ones(M_sectors*numUE,1); r_hi(~isIndoorUE) = params.deployRange;
for n = 1:nbrOfRealizations
    if n > 1
        %advance the trajectories by dt_snap; specular reflection keeps
        %each user inside its motion region (annulus around its site)
        UEpos_t = UEpos_t + v_move*params.dt_snap.*exp(1i*heading);
        r = abs(UEpos_t - siteCenterUE);
        out = (r > r_hi) | (r < r_lo);
        if any(out)
            u = exp(1i*heading(out));
            nvec = (UEpos_t(out) - siteCenterUE(out))./max(r(out),1e-9);
            u_ref = u - 2*real(u.*conj(nvec)).*nvec;
            heading(out) = angle(u_ref);
            r_cl = min(max(r(out), r_lo(out)), r_hi(out));
            UEpos_t(out) = siteCenterUE(out) + r_cl.*nvec;
        end
        params.UE_locations = [real(UEpos_t), imag(UEpos_t)];
        params.numUE = numUE;
        params.numCPE = numCPE_tot;
        params.CPE_locations = CPE_locations;
        %re-evaluate the setup at the new positions, evolving the latent
        %large-scale variables consistently (seed 0: the RNG stream
        %continues across snapshots)
        [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,~,D_FWA,D_cell,~,~,~,~,~,~,~,repDonorGaindB,mobState] = generateSetup(params,0,mobState);
        params.rep_donor_gain = db2pow(repDonorGaindB);
        params.gainOverNoise_lin = db2pow(gainOverNoisedB);
        params.BETA_interUE = db2pow(gainOverNoisedB_ue);
    end
    %(a) cellular UEs, NCR-assisted (proposed system)
    params.numUE = numUE;
    params.numCPE = numCPE_tot;
    params.CPE_locations = CPE_locations;
    params.D_FWA = D_FWA;
    params.D_cell = D_cell;
    params.R_gNB = R_gNB;
    params.R_cpe = R_cpe;
    params.R_interue = R_interue;
    params.R_ue = R_ue;
    params.set_repeat = set_repeat_pool;
    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
    rate_cell_rep(:,n) = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est);
    %(b) cellular UEs, plain OFDM (no NCRs)
    params.numCPE = 0;
    params.CPE_locations = [];
    params.D = D_cell;
    params.R_gNB = R_gNB(:,:,:,1+numCPE_tot:end);
    params.R_cpe = [];
    params.R_interue = [];
    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, ~, ~] = computePhysicalChannels_sub6_MIMO(params);
    rate_cell_plain(:,n) = compute_link_rates_OFDM(params, channel_dl, channel_est_dl);
    %(c) FWA CPEs (stationary: zero displacement freezes their setup)
    params.numUE = 0;
    params.numCPE = numCPE_tot;
    params.CPE_locations = CPE_locations;
    params.D_cell = [];
    params.D_FWA = D_FWA;
    params.D = D_FWA;
    params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
    params.R_cpe = R_cpe;
    params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
    params.R_ue = [];
    params.set_repeat = [];
    params.SI_cancel_factor = 10^(0.1*SI_cancel_dB);
    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, ~, ~] = computePhysicalChannels_sub6_MIMO(params);
    rate_FWA(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);
    %restore the full user set before the next snapshot's setup call
    params.numUE = numUE;
    params.R_ue = R_ue;
end

%% Save raw per-user x per-realization rate matrices for cross-seed pooling
dataFolder = fullfile(pwd,'resultData');
rateFolder = fullfile(dataFolder,'FWA_packing_analysis');
if not(isfolder(rateFolder))
    mkdir(rateFolder)
end
writematrix(rate_cell_rep,  fullfile(rateFolder,strcat('packing_cell_rep_',aID,'.csv')));
writematrix(rate_cell_plain,fullfile(rateFolder,strcat('packing_cell_plain_',aID,'.csv')));
writematrix(rate_FWA,       fullfile(rateFolder,strcat('packing_FWA_',aID,'.csv')));

%% In-run summary: safe load fraction f(eps) = q_eps(rate)/mean(rate).
%This seed only (temporal distribution); the cross-seed cellular pooling
%that adds population uncertainty happens in dataProcess/plot_packing.m
fprintf('\nseed %s: safe load fraction f = q_eps/mean (this drop only)\n', aID);
fprintf('%8s %12s %12s %12s %10s\n','eps','cell (NCR)','cell (plain)','FWA','multiple');
for e = eps_grid
    f_cr = quantile(rate_cell_rep(:),e)/mean(rate_cell_rep(:));
    f_cp = quantile(rate_cell_plain(:),e)/mean(rate_cell_plain(:));
    %FWA: admission is per-CPE with a measured link, so the usable
    %fraction is the ratio of summed per-CPE floors to summed means
    f_fw = sum(quantile(rate_FWA,e,2))/sum(mean(rate_FWA,2));
    fprintf('%8.2f %12.3f %12.3f %12.3f %10.2fx\n', e, f_cr, f_cp, f_fw, f_fw/f_cr);
end
toc(tStart)
