close all;
clear;
tStart = tic;
aID = getenv('SLURM_ARRAY_TASK_ID');

% This is for running on a cluster in parallel
% the bash script should give the aID as input
if (isempty(aID))
    warning('aID is empty. Trying SLURM ID.')
    aID = getenv('SLURM_ARRAY_TASK_ID');
end
if(isempty(aID))
    warning('aID is empty. Replacing it with 0010.')
    aID = '0084';
end
%RNG seed.
rng(str2double(aID),'twister');
%% Mobility params (TR 38.901 Table 7.2-5)
params.MOBILE = 1;
params.ue_velocity_outdoor = 40/3.6; %in-car UTs: 40 km/h (TR 38.901 Table 7.2-5, SMa)
params.ue_velocity_indoor = 3/3.6;   %indoor and pedestrian UTs: 3 km/h (Table 7.2-5)
%CSI aging lag: time between channel sounding and its use for precoding.
%C-band NR runs 30 kHz SCS (0.5 ms slots); with the typical 5 ms
%SRS/CSI-RS sounding periodicity (configurable per TS 38.214/TS 38.331)
%the precoder is on average ~2.5 ms stale when applied - this lag, not a
%single OFDM symbol, is what ages the channel for moving users
params.Ts = 2.5e-3; %CSI staleness (s): half the 5 ms sounding period
params.fc = 3.7e9; %carrier frequency (Hz)
params.c = 3e8; %speed of light (m/s)
%params.mob_rho (one channel-aging factor per cellular UE) is set after
%generateSetup, once the indoor/outdoor state of each UE is known
%% Rate requirements
%FCC 2024 Section 706 Report (FCC 24-27, https://docs.fcc.gov/public/attachments/FCC-24-27A1.pdf):
% - mobile 5G-NR coverage benchmark 35/3 Mbps DL/UL
% - fixed broadband benchmark 100/20 Mbps
% - long-term fixed broadband goal 1000/500 Mbps
params.r_min_cell = 35e6; %cellular DL floor at the 5th percentile: the FCC
                          %mobile 5G-NR coverage benchmark of 35/3 Mbps DL/UL
                          %(FCC 24-27); attainable here because the FWA/Wi-Fi
                          %offload carries the indoor demand
params.r_max_FWA = 1e9;   %FWA rate cap = FCC long-term fixed broadband goal (1 Gbps),
                          %also the Verizon 5G Home "Ultimate" gigabit tier
%%
params.BEAM = 0;
%%


%FWA CPE count per site, swept as a TAKE RATE: each site's apartment
%ring plus its share of the homes belt hold ~450 households, so 25-100
%CPEs spans a ~5-20%% take rate. US FWA serves ~12%% of
%broadband households in mid-2025, ~15%% by early 2026
%(https://www.telecompetitor.com/fixed-wireless-is-thriving-and-that-could-be-a-problem-report/;
%https://insights.opensignal.com/2025/10/20/the-state-of-us-fwa-what-impact-has-att-internet-airs-launch-had/dt;
%https://www.ericsson.com/en/reports-and-papers/mobility-report/dataforecasts/fwa-outlook)
numCPE_all_arr = 25:25:100; %CPEs per site (take-rate sweep)
params.cpe_mall_frac = 0.10; %fraction of all CPEs on the strip-mall rooftop
params.cpe_home_frac = 0.35; %fraction on single-family homes (FWA's core market)
params.cpe_op_frac = 0.20;   %fraction on office-park rooftops in the outer belt
                             %- the CPE hosts nearest the noise-limited edge;
                             %the remainder sit on apartment top floors
Band = 100e6;
params.noiseFigure = 9; %receiver noise figure (dB)

%% Suburban Macro (SMa) deployment, 3GPP TR 38.901 Table 7.2-5
%FWA's market case is strongest where fiber is scarce, i.e. suburbs -
%the urban (UMa) detour is reverted. Two macro sites are deployed
%CENTRALLY at (+/- ISD/2, 0) along a randomly oriented axis, so the
%deployment region surrounds them symmetrically.
params.ISD = 1732; %inter-site distance (m); the LARGER of the two SMa
                   %values listed in Table 7.2-5 (1299 / 1732 m), giving a
                   %coverage-limited cell edge where the NCR mechanisms can
                   %convert their link-budget gains into rate; the extra
                   %area fills with vegetation and roads (below)
params.deployRange = params.ISD/2; %cell radius (half ISD); sets the wrap-around span
%SUBURBAN ZONING, radially per site plus one shared midpoint zone:
% - apartment complex: ring [min_dist_2D, aptRadius] around each site,
%   2-4 floor buildings, FWA CPEs on top floors;
% - single-family homes: ring [aptRadius, homeRadius], 1-2 floors;
% - strip mall: disc of mallRadius at the ORIGIN - the arterial midpoint
%   between the two neighbourhoods, classic suburban commercial siting -
%   single-storey, high-loss commercial construction, FWA CPEs on the
%   rooftop (strategic NCR donors with open sightlines to both sites);
% - office parks: n_officepark small discs per site at opDist on the
%   outer roads, low-rise HIGH-LOSS commercial buildings (same M.2412
%   commercial-glass observation as the mall) hosting rooftop CPEs;
% - natural vegetation buffer: ring [homeRadius, deployRange], no users,
%   entering the SMa LOS model through the vegetation clutter density -
%   at ISD 1732 this belt spans 450-866 m, filling the larger cell;
% - roads: in-car users over [aptRadius, deployRange], whose area (and
%   hence user count) grows with the ISD.
params.aptRadius = 175;  %apartment-complex ring outer radius (m)
params.homeRadius = 450; %single-family belt outer radius (m)
params.mallRadius = 75;  %strip-mall disc radius (m), centred at the origin
%OFFICE PARKS: smaller commercial clusters on the arterial roads in the
%outer belt - low-rise office buildings with FWA CPEs on their rooftops,
%whose co-located NCRs sit close to the noise-limited road users at the
%coverage-limited cell edge (the CPE-hosted deployment logic of the
%paper: repeaters go where subscriber equipment already is)
params.n_officepark = 2; %office parks per site
params.opRadius = 60;    %office-park parcel radius (m)
params.opDist = 600;     %office-park centre distance from its site (m),
                         %on the roads between the homes belt and the edge
params.op_floors = [2 3]; %low-rise office buildings, floors
params.apt_floors = [2 4];  %apartment buildings, floors (suburban complexes)
params.home_floors = [1 2]; %single-family homes, floors
params.wrap_margin = 30; %guard (m) added to the wrap-around square so wrapped
                         %replicas never coincide with real positions
params.coverageRange_sub6 = 2*params.deployRange + params.wrap_margin; %wrap-around square length
params.min_dist_2D = 35; %minimum BS-UT 2D distance (m)
params.num_antennas_per_gNB = 64;
params.rho_tot = 10^(0.1*56); %BS Tx power 56 dBm per sector for the 100 MHz band:
                              %the TR 38.901 clause 7.8 SMa calibration assumptions
                              %specify "BS Tx power 49 dBm for SMa, bandwidth 20 MHz"
                              %(same PSD as the ITU-R M.2412 macro configurations);
                              %scaled at constant PSD to 100 MHz this gives 56 dBm,
                              %consistent
                              %with commercial 320 W C-band massive-MIMO AAUs. FCC
                              %compliance: peak-beam EIRP PSD = 36 dBm/MHz conducted
                              %+ 8 dBi element + 18.1 dB array gain (64 el.) =
                              %62.1 dBm/MHz, just under the 1640 W/MHz (62.15
                              %dBm/MHz) non-rural C-band limit of 47 CFR 27.50(j)
                              %(https://www.ecfr.gov/current/title-47/chapter-I/subchapter-B/part-27/subpart-C/section-27.50).
                              %CONSTANT-PSD CONVENTION: rho_tot is spread uniformly
                              %over the full band, so a service occupying a band
                              %fraction x radiates x*rho_tot and sees x of the
                              %noise - power splits in exact proportion to the
                              %band split between cellular and FWA, and every
                              %SINR computed with full-band power against
                              %full-band-normalized channels is exact on any
                              %subband (verified: halving band and power leaves
                              %all SINRs unchanged and halves all rates to
                              %machine precision)
params.CELL_REPEAT = 1;
params.TEMPORAL_MOBILITY = 1; %1 = cellular UEs move along trajectories and the
                              %large-scale state evolves across snapshots (see the
                              %mobility model below); 0 = one static setup per seed
%Number of antennas per UE
params.N_UE_FWA = 8;
params.N_UE_cell = 2; %4;
params.hr = 1.5; %outdoor UT height (m)
%CPE mounting height is drawn per CPE in generateSetup.m: top floor of
%its apartment building (params.apt_floors) or the strip-mall rooftop
params.ht_bs = 35; %BS antenna height (m), SMa (Table 7.2-5)
%% Sectored BS antennas: 3 sectors per site, element pattern per TR 38.901 Table 7.3-1
params.sectors_per_site = 3;
sector_offsets = [30; 150; 270]; %sector boresight azimuths within a site (deg)
params.theta_3dB = 65; %vertical 3 dB beamwidth (deg)
params.phi_3dB = 65; %horizontal 3 dB beamwidth (deg)
params.SLA_V = 30; %side-lobe attenuation, vertical (dB)
params.A_max = 30; %front-back attenuation (dB)
params.G_ant_max = 8; %max element gain (dBi)
%Mechanical downtilt in zenith-angle convention (90 deg = horizon):
%92 deg per the TR 38.901 clause 7.8 SMa calibration for ISD = 1732 m
%(95 deg for ISD = 1299 m)
params.tilt_zenith = 92;
%% SMa propagation, TR 38.901 Clause 7.4 (pathloss, LOS probability, O2I)
params.h_bldg = 10;  %average building height h (m), SMa default
params.W_street = 10; %average street width W (m), SMa default
params.r_vegetation = 0.20; %vegetation clutter density, set to DENSE - the
                            %highest of the three values defined in TR 38.901
                            %Table 7.4.2-1 (SMa row): 0 = no vegetation,
                            %0.10 = sparse, 0.20 = dense; suburban areas have
                            %"a high variability of foliage density" and more
                            %vegetation than urban (TR 38.901 clause 7.2)
params.indoor_UT_ratio = 0.8; %fraction of cellular UEs indoors (Table 7.2-5)
params.high_loss_ratio = 0; %HIGH-LOSS (IRR-glass) share of RESIDENTIAL
       %buildings: zero, because "the use of such high-loss glass
       %currently appears to be more predominant in commercial buildings
       %than in residential buildings" (ITU-R M.2412 / TR 38.901 Clause
       %7.4.3.1); the strip mall is commercial and always high-loss
       %(enforced by zone in generateSetup.m)
lambda_BS = 5;
lambda_SC = 0; %no small cells (kept for the result-file format)
%Active-user densities: GROSS busy-hour concurrently active data users
%on ANY access network, per suburban zone. Apartment complexes hold
%~4,000-6,000 residents/km^2 and single-family suburbs ~1,000-1,500
%(Demographia, http://demographia.com/db-intlsub.htm); a ~5% busy-hour
%active share - the top of the typical 3-5% suburban range, reflecting
%high-data-rate evening applications - gives ~250 and ~60 gross active
%users/km^2 respectively. The strip mall carries a modest evening crowd
%(~150/km^2 over its small parcel). Suburban roads carry ~4-15 active
%in-car users/km^2. The resulting ~12 cellular-active UEs per sector
%sits near the 10-per-TRxP full-buffer evaluation convention of ITU-R
%M.2412 / 3GPP TR 38.913.
%FWA/WI-FI OFFLOAD: cellular guarantees coverage to OUTDOOR users
%primarily - indoor users ride the home/venue connection (the FWA CPEs
%of this study). Wi-Fi/femtocell offload carries ~60% of ALL mobile data
%globally and 66% in the US (Cisco VNI,
%https://wballiance.com/the-cisco-visual-networking-index-vni-global-mobile-data-traffic-forecast-update-2017-2021/;
%https://newsroom.cisco.com/c/r/newsroom/en/us/a/y2019/m02/cisco-global-mobile-networks-will-support-more-than-12-billion-mobile-devices-and-iot-connections-by-2022-mobile-traffic-approaching-the-zettabyte-mil.html),
%concentrated at home where ~81% of US households have home internet
%(https://insights.opensignal.com/2025/10/20/the-state-of-us-fwa-what-impact-has-att-internet-airs-launch-had/dt);
%an 80% indoor offload share follows, and the offloaded indoor demand IS
%the FWA traffic this paper monetizes. Cellular user classes:
% - OUTDOOR PEDESTRIANS (3 km/h, no O2I, 1.5 m) in the apartment ring,
%   the homes belt, and the mall parking - a pedestrian_ratio share of
%   each zone's gross actives (TR 38.901 Table 7.2-5 indoor:outdoor);
% - IN-CAR users (40 km/h, Table 7.2-5) on the roads [aptRadius,
%   deployRange];
% - the 20% INDOOR remainder in each zone (3 km/h, zone-specific floor
%   heights and O2I class).
lambda_UE_apt_arr = 250; %125:125:500; %GROSS active users per km^2, apartment ring
lambda_UE_home = 60;     %GROSS active users per km^2, single-family belt
lambda_UE_mall = 150;    %GROSS active users per km^2, strip-mall parcel
lambda_UE_op = 200;      %GROSS active users per km^2, office-park parcels
                         %(sparse evening staff/visitors; the parks' main
                         %role in the busy hour is hosting CPEs/NCRs)
lambda_UE_road_arr = 10; %4:4:16;       %active in-car users per km^2, roads (paired)
params.pedestrian_ratio = 0.2; %fraction of each zone's GROSS actives outdoors on foot
params.indoor_offload = 0.8;   %fraction of indoor actives on home/venue broadband
                               %instead of cellular (Cisco VNI offload + US home
                               %internet penetration, see above)
params.Lmax = 1;
params.preLogFactor = 1;
params.loss_pc_cell = 5/100; %enforce r_min_cell at the 5th-percentile cell UE
%FWA multi-user/inter-cell interference suppression factor gamma_I (dB):
%residual interference power after beam nulling by the stationary,
%large-array FWA CPEs. Justified values: commercial ngFWA radios cancel
%interference by up to 45 dB (Tarana G1,
%https://taranawireless.com/principles-of-ngfwa-design-interference-cancellation/;
%ngFWA primer, https://resourcesapi.taranawireless.com/storage/resource_files/white-papers/1758649242_Tarana-ngFWA-Primer-White-Paper-2509-02.pdf)
%and real-time adaptive arrays measure ~30-35 dB null depths, still
%>20 dB with 3-6 simultaneous interferers (Sensors 2023,
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10386719/). The -20 dB
%default is conservative against both; sweep -20:5:0 for the
%suppression-factor SE comparison figure
SI_cancel_arr = -20; %-20:5:0;
params.HW_IMPAIRMENTS = 1;  % 1 = hardware impairments on, 0 = ideal hardware
%EVM-based impairment factors: an error vector magnitude of e turns a
%fraction e^2 of the signal power into distortion, so K = 1 - e^2.
%3GPP minimum signal quality: BS transmit EVM 3.5% for 256QAM / 8% for
%64QAM, TS 38.104 Table 6.5.2-1 (https://www.etsi.org/deliver/etsi_ts/138100_138199/138104/);
%UE EVM uses the same classes (TS 38.101-1); NR repeater conformance
%reuses them too (TS 38.115-2, https://www.etsi.org/deliver/etsi_ts/138100_138199/13811502/17.03.00_60/ts_13811502v170300p.pdf;
%RF requirements in TS 38.106). 256QAM-grade (3.5%) is assumed for the
%BS transmit chain, 64QAM-grade (8%) for the UE receive chain and for
%each half of the repeater amplify-and-forward chain.
params.Kt = 1 - 0.035^2;    % BS transmit chain, 256QAM EVM 3.5% (TS 38.104)
params.Kr = 1 - 0.08^2;     % UE receive chain, 64QAM-grade EVM 8%
params.Kt_rep = 1 - 0.08^2; % repeater transmit chain, EVM 8% (TS 38.115-2)
params.Kr_rep = 1 - 0.08^2; % repeater receive chain, EVM 8% (TS 38.115-2)
params.IMPERFECT_CSI = 1;   % 1 = beamforming uses MMSE channel estimates from noisy
                            %     CSI-RS observations, 0 = perfect CSI
params.csi_rs_offset_dB = 0;% CSI-RS EPRE relative to PDSCH EPRE (powerControlOffset,
                            % TS 38.214; standardised range -8..15 dB, default 0)
%Network-controlled repeater (NCR) amplification, per 3GPP: fixed
%amplification gain of 90 dB subject to the maximum output power
%constraint of 33 dBm (evaluation assumptions of TR 38.867, "Study on NR
%network-controlled repeaters"; RF repeater requirements in TS 38.106).
%The cap is applied per repeater at its actual operating point below, so
%the effective gain is min(G_max, P_max/P_in)
rep_gain_arr = 90; %NCR amplification gain G_max (dB), TR 38.867
params.rep_max_pow_dBm = 33; %NCR maximum output power (dBm), TR 38.867 / TS 38.106
params.ncr_zf_donor = 1; %1 = interference-aware (max-SINR / regularized
                         %zero-forcing) donor receive combining at the NCR,
                         %suppressing the co-channel sites at the donor input
                         %- the SIR ceiling of the repeated path; 0 = matched
                         %combining only (ablation)
params.ncr_rank = 1; %NCR forwarded spatial layers. 1 (default) = Rel-18 NCR:
                     %analog beamforming, a single amplify-and-forward chain,
                     %no digital processing (Carvalho et al., "Network-
                     %controlled repeater - an introduction", IEEE Comm.
                     %Standards Mag. 2026, https://arxiv.org/pdf/2403.09601).
                     %Set to params.N_UE_cell for the beyond-Rel-18
                     %"MIMO-capable" repeater what-if: one AF chain per UE
                     %antenna (top singular beam pairs, output power split
                     %across chains; clamped internally to the supported rank)
num_rep_arr = 1:1:6; %total enabled repeaters; greedy max-coverage attachment
                     %is per donor sector, saturating near one per sector (6)
params.num_repeater_per_cpe = 2; %NCRs attached per assisted user: TWO independent
                                 %rank-1 amplify-and-forward branches make the
                                 %composite channel rank-2, restoring the UE's
                                 %second spatial stream that a single Rel-18 NCR
                                 %collapses
params.rep_assist_frac = 0.2;    %fraction of each sector's WEAKEST cellular UEs (cell edge,
                                 %by serving-link large-scale gain) eligible for repeater
                                 %assistance; the rest keep the pure direct path
%Three layers of randomness, averaged per SLURM array task:
%  (1) the seed aID fixes the DROP: gNB/CPE/UE locations and the initial
%      large-scale state;
%  (2) nbrOfSnapshots units of MOBILITY: the cellular UEs advance
%      dt_snap along their trajectories and the large-scale state
%      (shadowing, LOS, association) evolves consistently. Motion model:
%      straight-line constant speed with a random heading per user
%      (Random Direction model: Camp, Boleng & Davies, WCMC 2002,
%      https://onlinelibrary.wiley.com/doi/10.1002/wcm.72), the 3GPP /
%      ITU-R M.2412 evaluation convention; in-car road UEs at 30 km/h in the
%      outer ring follow the straight-segment core of the Freeway
%      highway model (Bai, Sadagopan & Helmy, IEEE INFOCOM 2003,
%      https://ieeexplore.ieee.org/document/1208920), indoor UEs walk at
%      3 km/h inside the complex (TR 38.901 Table 7.2-1; in-car speed
%      per ITU-R M.2412 Dense Urban eMBB), with
%      specular reflection at the region boundaries;
%  (3) nbrOfRealizations fading / noise / CSI-error realizations per
%      snapshot (within a snapshot the CSI-lag Jakes channel aging
%      applies).
%Cellular statistics average over layers 2 x 3. The FWA CPEs are
%stationary (no mobility layer): they average over
%nbrOfSnapshots*nbrOfRealizations fading draws, so both services rest on
%the same sample count.
params.nbrOfSnapshots = 5;     %mobility units per seed (layer 2)
params.nbrOfRealizations = 10; %fading realizations per snapshot (layer 3)
params.dt_snap = 1;            %time between mobility snapshots (s)

%Minimum guaranteed FWA rate per subscriber-plan tier (swept). Values
%span the US 5G-FWA plan landscape (Q4 2024): 50 Mbps entry floor,
%100 Mbps = FCC fixed-broadband benchmark and common base plan,
%150-200 Mbps = current median delivered speeds (Verizon/AT&T ~150,
%T-Mobile ~205), 300 Mbps = Verizon 5G Home standard-plan soft cap;
%the gigabit "Ultimate" tier is captured by params.r_max_FWA = 1 Gbps.
%Refs: Light Reading (https://www.lightreading.com/fixed-wireless-access/verizon-and-t-mobile-fwa-speeds-are-actually-increasing),
%tecknexus (https://tecknexus.com/t-mobile-verizon-and-att-boost-5g-fwa-speeds-and-subscribers/),
%Verizon FWA guide (https://www.verizon.com/home/internet/guides/what-is-fixed-wireless-access-fwa-technology/).
r_min_arr = 1e6*(50:50:300);
%% Simulation FR1 setup
for idxBSDensity = 1:length(lambda_BS)
    %% gNB locations: macro sites, each split into 3 co-located sector BSs
    params.M_sites = 2; %number of macro sites. The TWO-SITE side-by-side
    %topology is structural, not just a parameter: the site placement
    %(Rsite = [0; ISD]), the one-complex-per-site layout, and the
    %wrap-around span all assume it - changing M_sites needs layout rework
    M_sites = params.M_sites;
    S = params.sectors_per_site;
    params.numGNB = M_sites*S; %each sector is an independent BS entry
    M_sectors = params.numGNB;
    %two sites centred about the origin at (+/- ISD/2) along a randomly
    %oriented axis; the strip mall sits at their midpoint (the origin)
    angleSite = 2*pi*rand;
    Rsite = [-params.ISD/2; params.ISD/2];
    siteLocations = [Rsite.*cos(angleSite), Rsite.*sin(angleSite)];
    params.siteLocations = siteLocations;
    params.locationsBS = kron(siteLocations, ones(S,1)); %sectors co-located at their site
    params.sector_boresights = repmat(sector_offsets, M_sites, 1); %boresight azimuth per BS entry (deg)
    for idxnumCPE = 1:length(numCPE_all_arr) %FWA take-rate sweep
    numCPE_all = numCPE_all_arr(idxnumCPE);
    %% CPE locations: numCPE_all per site across THREE zones - apartment
    %top floors, single-family rooftops in the homes belt (FWA's core
    %market), and the strip-mall rooftop at the origin
    numCPE_tot = M_sites*numCPE_all;
    %office-park centres: n_officepark per site at opDist, random azimuths
    opC = zeros(M_sites*params.n_officepark,1);
    for idxSite = 1:M_sites
        thOP = 2*pi*rand(params.n_officepark,1);
        opC((idxSite-1)*params.n_officepark+(1:params.n_officepark)) = ...
            siteLocations(idxSite,1) + 1i*siteLocations(idxSite,2) + params.opDist*exp(1i*thOP);
    end
    n_parks = numel(opC);
    n_cpe_mall = round(params.cpe_mall_frac*numCPE_tot);
    n_cpe_op = round(params.cpe_op_frac*numCPE_tot);
    n_home_site = round(params.cpe_home_frac*numCPE_tot/M_sites);
    n_apt_site = round((numCPE_tot - n_cpe_mall - n_cpe_op - M_sites*n_home_site)/M_sites);
    n_cpe_op = numCPE_tot - M_sites*(n_apt_site + n_home_site) - n_cpe_mall;
    RCPE = sqrt(params.min_dist_2D^2 + (params.aptRadius^2 - params.min_dist_2D^2)*rand(M_sites*n_apt_site,1));
    angleCPE = 2*pi*rand(M_sites*n_apt_site,1);
    CPE_apt = kron(siteLocations, ones(n_apt_site,1)) + [RCPE.*cos(angleCPE), RCPE.*sin(angleCPE)];
    RHme = sqrt(params.aptRadius^2 + (params.homeRadius^2 - params.aptRadius^2)*rand(M_sites*n_home_site,1));
    angleHme = 2*pi*rand(M_sites*n_home_site,1);
    CPE_home = kron(siteLocations, ones(n_home_site,1)) + [RHme.*cos(angleHme), RHme.*sin(angleHme)];
    RMall = params.mallRadius*sqrt(rand(n_cpe_mall,1));
    angleMall = 2*pi*rand(n_cpe_mall,1);
    CPE_mall = [RMall.*cos(angleMall), RMall.*sin(angleMall)]; %at the origin
    %office-park CPEs: round-robin across the parks, uniform in each parcel
    parkOfCPE = 1 + mod(0:n_cpe_op-1, n_parks)';
    posOP = opC(parkOfCPE) + params.opRadius*sqrt(rand(n_cpe_op,1)).*exp(2i*pi*rand(n_cpe_op,1));
    CPE_op = [real(posOP), imag(posOP)];
    CPE_locations = [CPE_apt; CPE_home; CPE_mall; CPE_op];
    %CPE zone ids: 1 apartment, 2 single-family home, 3 strip mall, 4 office park
    params.cpe_zone = [ones(M_sites*n_apt_site,1); 2*ones(M_sites*n_home_site,1); ...
                       3*ones(n_cpe_mall,1); 4*ones(n_cpe_op,1)];
    for idxUEDensity = 1:length(lambda_UE_apt_arr)
        lambda_road = lambda_UE_road_arr(min(idxUEDensity,end));
        %% UE construction, per suburban zone. For each zone the GROSS
        %actives split into outdoor pedestrians (pedestrian_ratio, all
        %cellular) and indoor users, of which only (1 - indoor_offload)
        %remain cellular. Every cellular UE carries its zone id, motion
        %region (centre + annulus) and speed, built here explicitly.
        A_apt = pi*((params.aptRadius/1000)^2 - (params.min_dist_2D/1000)^2);  %km^2 per site
        A_home = pi*((params.homeRadius/1000)^2 - (params.aptRadius/1000)^2);  %km^2 per site
        A_mall = pi*(params.mallRadius/1000)^2;                                %km^2, shared
        A_roads = pi*((params.deployRange/1000)^2 - (params.aptRadius/1000)^2);%km^2 per site
        g_apt = ceil(lambda_UE_apt_arr(idxUEDensity)*A_apt);   %gross actives per site
        g_home = ceil(lambda_UE_home*A_home);                  %gross actives per site
        g_mall = ceil(lambda_UE_mall*A_mall);                  %gross actives, shared zone
        A_park = pi*(params.opRadius/1000)^2;                  %km^2 per office park
        g_park = ceil(lambda_UE_op*A_park);                    %gross actives per park
        n_car = ceil(lambda_road*A_roads);                     %in-car per site (all cellular)
        ped_apt = round(params.pedestrian_ratio*g_apt);
        ped_home = round(params.pedestrian_ratio*g_home);
        ped_mall = round(params.pedestrian_ratio*g_mall);
        ped_park = round(params.pedestrian_ratio*g_park);
        ind_apt = round((1 - params.indoor_offload)*(g_apt - ped_apt));
        ind_home = round((1 - params.indoor_offload)*(g_home - ped_home));
        ind_mall = round((1 - params.indoor_offload)*(g_mall - ped_mall));
        ind_park = round((1 - params.indoor_offload)*(g_park - ped_park));
        tot_cell = M_sites*(ped_apt + ind_apt + ped_home + ind_home + n_car) ...
                   + ped_mall + ind_mall + n_parks*(ped_park + ind_park);
        %total padded UP to a multiple of the sector count with extra
        %single-family pedestrians (the largest, least sensitive class)
        params.numUE = ceil(tot_cell/M_sectors);
        numUE = params.numUE;
        n_pad = M_sectors*numUE - tot_cell;
        %per-UE records: position, zone (1 apt / 2 home / 3 mall),
        %indoor flag, motion-region centre and annulus, speed
        UE_locations = zeros(M_sectors*numUE,2);
        ue_zone = ones(M_sectors*numUE,1);
        isIndoorUE = false(M_sectors*numUE,1);
        regionCenterUE = zeros(M_sectors*numUE,1); %complex, filled below
        r_lo_ue = zeros(M_sectors*numUE,1);
        r_hi_ue = zeros(M_sectors*numUE,1);
        v_move_ue = params.ue_velocity_indoor*ones(M_sectors*numUE,1);
        sitesC = siteLocations(:,1) + 1i*siteLocations(:,2);
        idx = 0;
        dropRing = @(c,rin,rout,n) c + sqrt(rin^2 + (rout^2 - rin^2)*rand(n,1)).*exp(2i*pi*rand(n,1));
        for idxSite = 1:M_sites
            zdef = { %count, zone, indoor, centre, r_lo, r_hi, speed
                ped_apt, 1, false, sitesC(idxSite), params.min_dist_2D, params.aptRadius, params.ue_velocity_indoor;
                ind_apt, 1, true,  sitesC(idxSite), params.min_dist_2D, params.aptRadius, params.ue_velocity_indoor;
                ped_home,2, false, sitesC(idxSite), params.aptRadius,   params.homeRadius,params.ue_velocity_indoor;
                ind_home,2, true,  sitesC(idxSite), params.aptRadius,   params.homeRadius,params.ue_velocity_indoor;
                n_car,   2, false, sitesC(idxSite), params.aptRadius,   params.deployRange,params.ue_velocity_outdoor};
            for z = 1:size(zdef,1)
                n_z = zdef{z,1};
                if n_z == 0, continue, end
                pos = dropRing(zdef{z,4}, zdef{z,5}, zdef{z,6}, n_z);
                UE_locations(idx+(1:n_z),:) = [real(pos), imag(pos)];
                ue_zone(idx+(1:n_z)) = zdef{z,2};
                isIndoorUE(idx+(1:n_z)) = zdef{z,3};
                regionCenterUE(idx+(1:n_z)) = zdef{z,4};
                r_lo_ue(idx+(1:n_z)) = zdef{z,5};
                r_hi_ue(idx+(1:n_z)) = zdef{z,6};
                v_move_ue(idx+(1:n_z)) = zdef{z,7};
                idx = idx + n_z;
            end
        end
        %shared strip-mall zone at the origin (5 m inner radius keeps the
        %reflection normal well-defined), then the single-family padding
        zdef = {
            ped_mall,3, false, 0, 5, params.mallRadius, params.ue_velocity_indoor;
            ind_mall,3, true,  0, 5, params.mallRadius, params.ue_velocity_indoor;
            n_pad,   2, false, sitesC(1+mod(0:n_pad-1,M_sites)).', params.aptRadius, params.homeRadius, params.ue_velocity_indoor};
        %office-park users, one block per park (zone 4)
        for pk = 1:n_parks
            zdef = [zdef;
                {ped_park,4, false, opC(pk), 5, params.opRadius, params.ue_velocity_indoor;
                 ind_park,4, true,  opC(pk), 5, params.opRadius, params.ue_velocity_indoor}]; %#ok<AGROW>
        end
        for z = 1:size(zdef,1)
            n_z = zdef{z,1};
            if n_z == 0, continue, end
            ctr = zdef{z,4};
            if isscalar(ctr), ctr = ctr*ones(n_z,1); else, ctr = ctr(:); end
            pos = ctr + sqrt(zdef{z,5}^2 + (zdef{z,6}^2 - zdef{z,5}^2)*rand(n_z,1)).*exp(2i*pi*rand(n_z,1));
            UE_locations(idx+(1:n_z),:) = [real(pos), imag(pos)];
            ue_zone(idx+(1:n_z)) = zdef{z,2};
            isIndoorUE(idx+(1:n_z)) = zdef{z,3};
            regionCenterUE(idx+(1:n_z)) = ctr;
            r_lo_ue(idx+(1:n_z)) = zdef{z,5};
            r_hi_ue(idx+(1:n_z)) = zdef{z,6};
            v_move_ue(idx+(1:n_z)) = zdef{z,7};
            idx = idx + n_z;
        end
        params.UE_locations = UE_locations;
        params.isIndoorUE = isIndoorUE;
        params.ue_zone = ue_zone;
        params.numCPE = numCPE_tot;
        params.CPE_locations = CPE_locations;
        params.Band = Band; %Communication bandwidth
        [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,pilotIndex,D_FWA,D_cell,APpositions,UEpositions,distances,distancesCPEs,isIndoor,hUT,O2IdB,repDonorGaindB,mobState] = generateSetup(params,str2double(aID));
        params.rep_donor_gain = db2pow(repDonorGaindB); %donor-side repeater sector antenna gain (linear power, M_sectors x K_FWA)
        params.gainOverNoise_lin = db2pow(gainOverNoisedB); %BS-to-user large-scale gain (linear), single source of truth
        %Total downlink power arriving at each repeater's donor panel (mW),
        %summed over all sector BSs: sets the operating point of the NCR
        %output-power cap. gainOverNoise_lin is normalized by the full-band
        %noise power, which is re-applied here
        noise_mW = 10^(0.1*(-174 + 10*log10(Band) + params.noiseFigure));
        P_in_rep_mW = params.rho_tot*noise_mW*sum(params.gainOverNoise_lin(:,1:numCPE_tot).*params.rep_donor_gain,1)';
        %Per-UE channel aging factor from the per-UE speeds built above
        %(pedestrians and indoor users 3 km/h, in-car 40 km/h)
        params.mob_rho = besselj(0,2*pi*params.Ts*v_move_ue*params.fc/params.c);
        params.BETA_interUE = db2pow(gainOverNoisedB_ue);
        %% Enabled-repeater pool order: GREEDY MAXIMUM COVERAGE over the
        %DISADVANTAGED cellular UEs. Eligibility mirrors the rate
        %computation exactly (the weakest rep_assist_frac fraction of
        %each sector's UEs by serving-link large-scale gain); a repeater
        %can cover an eligible UE iff its donor sector is the UE's
        %serving sector and the donor x service link metric is positive.
        %Each greedy pick covers the most NOT-YET-COVERED eligible UEs
        %(ties broken by the larger summed metric), and the pool for any
        %budget is a prefix of the greedy sequence. Greedy maximization
        %of a coverage (submodular) function is within (1 - 1/e) of the
        %optimal selection: G. L. Nemhauser, L. A. Wolsey, M. L. Fisher,
        %"An analysis of approximations for maximizing submodular set
        %functions - I", Mathematical Programming 14, 265-294 (1978),
        %https://doi.org/10.1007/BF01588971
        BETA_FWA_assoc = db2pow(gainOverNoisedB(:,1:numCPE_tot)).*D_FWA;
        donor_of_cpe = zeros(numCPE_tot,1);
        for q = 1:numCPE_tot
            donor_of_cpe(q) = find(D_FWA(:,q)==1,1);
        end
        eligible = false(M_sectors*numUE,1);
        gain_ue_drop = params.gainOverNoise_lin(:,numCPE_tot+1:end);
        for m = 1:M_sectors
            served = find(D_cell(m,:)==1);
            [~,ord] = sort(gain_ue_drop(m,served),'ascend');
            n_assist = ceil(params.rep_assist_frac*numel(served));
            eligible(served(ord(1:n_assist))) = true;
        end
        elig_idx = find(eligible)';
        Vcov = zeros(numCPE_tot,numel(elig_idx)); %coverage metric matrix
        for j = 1:numel(elig_idx)
            kUE = elig_idx(j);
            m_serv = find(D_cell(:,kUE)==1,1);
            v_rep = BETA_FWA_assoc(m_serv,:)'.*params.rep_donor_gain(m_serv,:)'.*params.BETA_interUE(1:numCPE_tot,numCPE_tot+kUE);
            v_rep(donor_of_cpe ~= m_serv) = 0;
            Vcov(:,j) = max(v_rep,0);
        end
        covered = false(1,numel(elig_idx));
        remaining = true(numCPE_tot,1);
        greedy_order = zeros(numCPE_tot,1);
        for pick = 1:numCPE_tot
            newcov = sum(Vcov(:,~covered) > 0, 2);
            newcov(~remaining) = -1;
            cand = find(newcov == max(newcov) & remaining);
            [~,ic] = max(sum(Vcov(cand,:),2)); %tie-break: summed metric
            greedy_order(pick) = cand(ic);
            covered = covered | (Vcov(cand(ic),:) > 0);
            remaining(cand(ic)) = false;
        end
        %% Trajectory of the cellular UEs (randomness layer 2)
        %The per-snapshot large-scale states are generated ONCE per seed,
        %before the repeater sweep, so every num_rep value sees the same
        %motion. Only the small gain/association matrices are stored; the
        %spatial correlation tensors (scaled identities) are rebuilt per
        %snapshot in the cell phase. The FWA CPEs have zero displacement,
        %so their rows/columns are frozen exactly by the evolution.
        nbrOfSnapshots = params.nbrOfSnapshots;
        nbrOfRealizations = params.nbrOfRealizations;
        snapGain = cell(nbrOfSnapshots,1);
        snapGainUE = cell(nbrOfSnapshots,1);
        snapDcell = cell(nbrOfSnapshots,1);
        %per-snapshot spatial-correlation CHOLESKY FACTORS: under the
        %local-scattering model they depend on the users' angles, so
        %they are captured from generateSetup per snapshot (R_cpe and
        %R_ue are frozen: CPEs are static and handsets are uncorrelated).
        %Memory: ~0.1-0.5 GB across snapshots at the largest CPE counts.
        snapR_gNB = cell(nbrOfSnapshots,1);
        snapR_int = cell(nbrOfSnapshots,1);
        snapGain{1} = gainOverNoisedB;
        snapGainUE{1} = gainOverNoisedB_ue;
        snapDcell{1} = D_cell;
        snapR_gNB{1} = R_gNB;
        snapR_int{1} = R_interue;
        if params.TEMPORAL_MOBILITY && nbrOfSnapshots > 1
            UEpos_t = UE_locations(:,1) + 1i*UE_locations(:,2);
            %per-UE motion regions and speeds were built with the zones
            siteCenterUE = regionCenterUE;
            v_move = v_move_ue;
            heading = 2*pi*rand(M_sectors*numUE,1);
            r_lo = r_lo_ue;
            r_hi = r_hi_ue;
            for t_snap = 2:nbrOfSnapshots
                %advance by dt_snap; specular reflection keeps each user
                %inside its motion region (annulus around its site)
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
                %seed 0: the RNG stream continues across snapshots
                [snapGain{t_snap},snapGainUE{t_snap},snapR_gNB{t_snap},~,snapR_int{t_snap},~,~,~,snapDcell{t_snap},~,~,~,~,~,~,~,~,mobState] = generateSetup(params,0,mobState);
            end
            params.UE_locations = UE_locations; %restore the drop positions
        else
            %static mode: every snapshot reuses the drop's setup
            for t_snap = 2:nbrOfSnapshots
                snapGain{t_snap} = gainOverNoisedB;
                snapGainUE{t_snap} = gainOverNoisedB_ue;
                snapDcell{t_snap} = D_cell;
                snapR_gNB{t_snap} = R_gNB;
                snapR_int{t_snap} = R_interue;
            end
        end
        for idxnumrep = 1:length(num_rep_arr)  
            for idxrepgain = 1:length(rep_gain_arr)
                params.repeat_gain = rep_gain_arr(idxrepgain); %G_max (dB), recorded in the CSV
                %full-band donor input power per repeater; the rate
                %function derives the SUBBAND-SELECTIVE effective gain
                %from it after attachment (TR 38.867 per-resource
                %side control)
                params.rep_P_in_mW = P_in_rep_mW;
                params.num_repeater_tot = num_rep_arr(idxnumrep); %total repeaters enabled
                %pool = prefix of the greedy maximum-coverage sequence
                %computed once per drop above
                params.set_repeat = greedy_order(1:min(params.num_repeater_tot,numCPE_tot));
                Band = params.Band; %Communication bandwidth
                K_FWA = params.numCPE;
                params.numUE = numUE;
                K = params.numCPE + M_sectors*params.numUE; 
                %Cell phase: mobility snapshots (layer 2) x fading
                %realizations (layer 3); the per-snapshot correlation
                %tensors are rebuilt from the stored gains
                rate_dl = zeros(K-K_FWA,nbrOfSnapshots*nbrOfRealizations);
                if params.CELL_REPEAT
                    params.D_FWA = D_FWA;
                    params.R_cpe = R_cpe;
                    params.R_ue = R_ue; 
                    for t_snap = 1:nbrOfSnapshots
                        params.gainOverNoise_lin = db2pow(snapGain{t_snap});
                        params.BETA_interUE = db2pow(snapGainUE{t_snap});
                        params.D_cell = snapDcell{t_snap};
                        params.R_gNB = snapR_gNB{t_snap};
                        params.R_interue = snapR_int{t_snap};
                        for n = 1:nbrOfRealizations
                            [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                            rate_dl(:,(t_snap-1)*nbrOfRealizations+n) = compute_link_rates_OFDM_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est);         
                        end
                    end
                else
                    params.numCPE = 0;
                    params.CPE_locations = [];
                    params.R_cpe = [];
                    params.R_interue = [];
                    params.R_ue = R_ue; 
                    for t_snap = 1:nbrOfSnapshots
                        params.gainOverNoise_lin = db2pow(snapGain{t_snap});
                        params.BETA_interUE = db2pow(snapGainUE{t_snap});
                        params.D = snapDcell{t_snap};
                        params.R_gNB = snapR_gNB{t_snap}(:,:,:,1+numCPE_tot:end);
                        for n = 1:nbrOfRealizations
                            [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                            rate_dl(:,(t_snap-1)*nbrOfRealizations+n) = compute_link_rates_OFDM(params, channel_dl, channel_est_dl);         
                        end
                    end
                end
                %restore the snapshot-1 (drop) state for the FWA phase; the
                %CPE rows/columns are identical in every snapshot anyway
                params.gainOverNoise_lin = db2pow(snapGain{1});
                params.BETA_interUE = db2pow(snapGainUE{1});
                mean_rate_dl_cell = mean(rate_dl,2);
                %Raw per-user x per-(snapshot x realization) cellular rates
                %for the packing-analysis figure (dataProcess/combineData.m);
                %the safe-load fraction f = q_eps/mean is bandwidth-scale-
                %invariant, so these full-band rates suffice
                packFolder = 'resultData/FWA_packing_analysis';
                if not(isfolder(packFolder))
                    mkdir(packFolder)
                end
                if params.CELL_REPEAT
                    cellTag = 'packing_cell_rep_';
                else
                    cellTag = 'packing_cell_plain_';
                end
                writematrix(rate_dl, strcat(packFolder,'/',cellTag,num2str(numCPE_all),'CPE_',num2str(params.num_repeater_tot),'rep_',aID,'.csv'));
                if (quantile(mean_rate_dl_cell,params.loss_pc_cell)>=params.r_min_cell)
                    params.Band = max(params.Band*(1-params.r_min_cell/quantile(mean_rate_dl_cell,params.loss_pc_cell)),0);
                    if (params.Band > 0)
                        mean_rate_dl_cell = mean_rate_dl_cell*params.r_min_cell/quantile(mean_rate_dl_cell,params.loss_pc_cell); 
                    end
                else
                    params.Band = 0;
                end
                Band_after_cell = params.Band; %band left after serving the cellular UEs
                %gamma_I sweep for the FWA phase (single -20 dB value by
                %default; see SI_cancel_arr for the sweep and citations)
                for idxSI = 1:length(SI_cancel_arr)
                    params.SI_cancel_factor = 10^(0.1*SI_cancel_arr(idxSI));
                    params.Band = Band_after_cell;
                    params.numUE = 0;
                    params.numCPE = numCPE_tot;
                    K_FWA = params.numCPE;
                    K = params.numCPE + M_sectors*params.numUE;
                    params.CPE_locations = CPE_locations;
                    params.D_cell = [];
                    params.D_FWA = D_FWA;
                    params.D = D_FWA;
                    CPE_idxs = 1:numCPE_tot;
                    params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
                    params.R_cpe = R_cpe;
                    params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
                    params.R_ue = []; 
                    params.set_repeat = [];
                    %stationary CPEs: no mobility layer; the fading draw
                    %count matches the cellular sample count (layers 2 x 3)
                    nbrFWADraws = nbrOfSnapshots*nbrOfRealizations;
                    rate_dl = zeros(K,nbrFWADraws);
                    for n = 1:nbrFWADraws
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    save_old_rate = rate_dl;
                    save_old_mean_FWA = mean_rate_dl_FWA;
                    if idxSI == 1 && params.Band > 0
                        %raw FWA packing matrix (skipped if the cell phase
                        %consumed the whole band: all-zero rates carry no
                        %distribution information)
                        writematrix(rate_dl, strcat(packFolder,'/packing_FWA_',num2str(numCPE_all),'CPE_',num2str(params.num_repeater_tot),'rep_',aID,'.csv'));
                    end
                    Band_FWA = params.Band;
                    for idxrmin = 1:length(r_min_arr)
                        params.r_min_FWA = r_min_arr(idxrmin);
                        K_FWA_max = 0;
                        params.Band = Band_FWA;
                        [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                        K_FWA_init = sum(FWA_util>0);
                        if (params.Band > 0)
                            while any(mean_rate_dl_FWA > params.r_max_FWA)
                                CPE_idxs = find(mean_rate_dl_FWA > params.r_max_FWA);
                                params.Band = params.Band*(1-(params.r_min_FWA/min(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA))));
                                rate_dl(CPE_idxs,:) = rate_dl(CPE_idxs,:)*params.r_min_FWA/min(mean_rate_dl_FWA(mean_rate_dl_FWA > params.r_max_FWA));
                                mean_rate_dl_FWA(CPE_idxs) = mean(rate_dl(CPE_idxs,:),2);
                                params.set_repeat = [params.set_repeat; CPE_idxs];
                                not_set_repeat = setdiff(1:numCPE_tot,params.set_repeat);
                                for n = 1:nbrFWADraws
                                    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, ~, ~] = computePhysicalChannels_sub6_MIMO(params);
                                    rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);                                              
                                end
                                mean_rate_dl_FWA(not_set_repeat) = mean(rate_dl(not_set_repeat,:),2);
                            end
                            [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                            K_FWA_max = sum(FWA_util>0);
                        end
                        sum_FWA_rate = sum(mean_rate_dl_FWA);
                        params.CPE_locations = CPE_locations;
                        params.numCPE = numCPE_tot;
                        params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
                        params.R_cpe = R_cpe;
                        params.R_interue = R_interue;
                        params.R_ue = []; 
                        params.D_FWA = D_FWA;
                        params.set_repeat = [];
                        params.Band = Band;
                        mean_rate_dl_FWA = save_old_mean_FWA;
                        rate_dl = save_old_rate;
                        if (Band_FWA > 0) && (K_FWA_max > 0)
                            FWA_se_out = sum_FWA_rate/(Band_FWA*K_FWA_max);
                        else
                            FWA_se_out = 0; %no band or no served CPEs: report 0 instead of 0/0
                        end

                    %% Recording the Results
                    
                        %Taking care of folder directory creation etc
                        dataFolder = 'resultData';
                        rateFolder = strcat(dataFolder,'/FWA_multi_cell_repeater_fix_comp_alloc');
                        if not(isfolder(dataFolder))
                            mkdir(dataFolder)
                        end
                        if not(isfolder(rateFolder))
                            mkdir(rateFolder)
                        end
                
                        deployRange = params.deployRange;
                        numBS = size(params.locationsBS,1);
                        result_string = strcat('/results_numFWA_',num2str(params.Lmax), 'Lmax_', num2str(numCPE_all),...
                            'CPE_',num2str(lambda_BS(idxBSDensity)),...
                            'lambdaBS_',num2str(lambda_SC(idxBSDensity)),...
                            'lambdaSC_',num2str(lambda_UE_apt_arr(idxUEDensity)),...
                            'lambdaUEres_',num2str(lambda_road),'lambdaUEroad_', ...
                            num2str(deployRange),'deployRange_', ...
                            num2str(params.r_min_FWA/10^6),'min_FWA_rate', ...
                            num2str(params.repeat_gain),'repeater_gain', ...
                            num2str(params.num_repeater_tot),'num_repeater_', ...
                            num2str(SI_cancel_arr(idxSI)),'SI_cancel_', ...
                            aID);
                        recording_text_file_string = strcat(rateFolder,result_string,'.csv');
                        fileID = fopen(recording_text_file_string,'w');
                        output_categories = ['lambdaBS,','lambdaSC,','numCPE,','lambdaUE,',...
                        'lambdaUE_road,','deployRange,','r_min_cell,','r_min_FWA,','num_rep,','rep_gain,','SI_cancel,','init_FWA,','max_FWA,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
                        fprintf(fileID,output_categories);
                        formatSpec = '%d,%d,%d,%d,%d,%d,%f,%f,%d,%d,%d,%d,%d,%f,%f,%f,%f\n';
                        fprintf(fileID,formatSpec,lambda_BS(idxBSDensity),lambda_SC(idxBSDensity),numCPE_all, ...
                        lambda_UE_apt_arr(idxUEDensity),lambda_road,deployRange,params.r_min_cell,params.r_min_FWA,params.num_repeater_tot,params.repeat_gain,SI_cancel_arr(idxSI),K_FWA_init,K_FWA_max,Band,Band_FWA, sum(mean_rate_dl_cell)/(Band - Band_FWA), FWA_se_out);
                        fclose(fileID);
                        % 'lambdaUE_road,','deployRange,','r_min_cell,','r_min_FWA,','num_rep,','rep_gain,','init_FWA,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n']; %'max_pow_fac,'','max_cell_util,','max_FWA_util
                        % lambda_UE_apt_arr(idxUEDensity),lambda_road,deployRange,params.r_min_cell,params.r_min_FWA,params.num_repeater_tot,params.repeat_gain,K_FWA_init,Band,Band_FWA, sum(mean_rate_dl_cell)/(Band - Band_FWA), sum_FWA_rate/(Band_FWA*K_FWA_init));
                    end
                end %idxSI (SI-cancellation sweep)
            end
        end
    end
    end %idxnumCPE (FWA take-rate sweep)
end
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)

