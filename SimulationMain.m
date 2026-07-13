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
params.enable_fwa_cap_realloc = 0; %gate for the r_max cap-redistribution loop in the
                          %FWA phase (K_FWA_max). DISABLED pending review: with
                          %demands near the cap, one exceeder's max(rmin/rate)
                          %can claim up to half the shared band per iteration,
                          %making K_FWA_max cliff-sensitive to ~1% rate changes
%%
params.BEAM = 0;
%%


%FWA CPE deployment: potential CPE SITES follow the BUILDING INVENTORY
%(one per apartment block, single-family home, office building, and mall
%unit), and the sweep is the TAKE RATE applied to that inventory. US FWA
%serves ~12%% of broadband households in mid-2025, ~15%% by early 2026
%(https://www.telecompetitor.com/fixed-wireless-is-thriving-and-that-could-be-a-problem-report/;
%https://insights.opensignal.com/2025/10/20/the-state-of-us-fwa-what-impact-has-att-internet-airs-launch-had/dt;
%https://www.ericsson.com/en/reports-and-papers/mobility-report/dataforecasts/fwa-outlook),
%so 5-20%% brackets the market. The inventory derives from the RESIDENT
%densities and household size below, keeping CPE supply proportional to
%the population it serves (previously the CPE count was static while the
%user population grew - the inconsistency this fixes).
take_rate_arr = 0.05:0.05:0.20; %RESIDENTIAL FWA take rate applied to the
                                %apartment-block and single-family inventory
params.take_rate_biz = 0.95; %BUSINESS take rate for the commercial inventory
       %(office buildings, mall units): "nearly all businesses (95%)
       %report having at least one broadband Internet connection" (FCC
       %Business Broadband Capability Survey,
       %https://www.fcc.gov/document/business-broadband-capability-survey-results;
       %see also https://advocacy.sba.gov/2026/01/13/issue-brief-no-23-small-businesses-access-to-broadband-internet/);
       %in this FWA-served deployment those connections are the CPEs
params.hh_size = 2.5;  %persons per household (US Census ~2.5)
params.n_apt_blocks = 24; %apartment blocks per complex (one potential
                          %building CPE each, serving the block as an MDU)
params.n_op_bldg = 3;  %office buildings per office park (one potential CPE each)
params.n_mall_units = 6; %strip-mall anchor units (one potential CPE each, shared zone)
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
%Two FIXED macro sites (M_sites below): the legacy BS-density sweep and
%its lambdaBS/lambdaSC result-file tags are removed
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
%resident densities (the demographic source for BOTH actives and the
%CPE building inventory); the swept quantity below is the busy-hour
%DEVICE ACTIVITY share, not the residents
params.rho_res_apt = 5000;  %residents per km^2, apartment ring
params.rho_res_home = 1200; %residents per km^2, single-family belt
params.activity_ref = 0.05; %REFERENCE busy-hour activity share: the per-zone
                         %gross-active densities below are the values observed
                         %at this share (apartment ring = rho_res_apt x 5%)
activity_arr = 0.025:0.025:0.1; %SWEPT busy-hour device-activity share; every
                         %zone's gross-active density (and the road density)
                         %scales linearly with activity_arr/activity_ref
lambda_UE_apt_base = 250;  %GROSS active users per km^2 at activity_ref, apartment ring
lambda_UE_home_base = 60;  %GROSS active users per km^2 at activity_ref, single-family belt
lambda_UE_mall_base = 150; %GROSS active users per km^2 at activity_ref, strip-mall parcel
lambda_UE_op_base = 200;   %GROSS active users per km^2 at activity_ref, office-park
                         %parcels (sparse evening staff/visitors; the parks' main
                         %role in the busy hour is hosting CPEs/NCRs)
lambda_UE_road_base = 10;  %active in-car users per km^2 at activity_ref, roads
params.pedestrian_ratio = 0.2; %fraction of each zone's GROSS actives outdoors on foot
params.indoor_offload = 0.8;   %fraction of indoor actives on home/venue broadband
                               %instead of cellular (Cisco VNI offload + US home
                               %internet penetration, see above)
%params.Lmax = 1; %RETIRED: the sectored deployment always serves each
                  %user from its single strongest sector; the legacy
                  %cell-free Lmax multi-AP association was removed from
                  %generateSetup.m
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
%value is conservative against both (a -20:5:0 sweep produced the legacy
%suppression-factor SE-comparison figure; the factor is FIXED now)
params.SI_cancel_dB = -20;
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
params.repeat_gain = 90; %FIXED NCR amplification gain G_max (dB), TR 38.867 (not swept)
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
params.ncr_benefit_gate = 1; %per-UE NCR ON/OFF side control (TR 38.867); =0 only for ablation
                     %is per donor sector, saturating near one per sector (6)
params.FWA_REPEAT = 1; %1 = HYBRID FWA scheduling with NCR assistance:
params.rep_assist_frac_fwa = 1; %FWA-side NCR admission: ALL needy CPEs are
                       %admitted for NCR support, DECOUPLED from the swept
                       %cellular rep_assist_frac (which only shapes cellular
                       %eligibility and the greedy pool)
                       %CPEs that miss their demand on the shared full band
                       %("needy") get a DEDICATED subband - the satisfied
                       %CPEs' shared subband shrinks just enough that every
                       %one still meets its demand (the binding CPE sits
                       %exactly at its demand), and the freed remainder
                       %serves the needy CPEs (who ALSO keep their f_sat
                       %share of the shared band) with NCR-aided MU-MIMO
                       %(compute_link_rates_MIMO_mmse_wi_repeater.m: same
                       %pool, donor-sector lock and benefit gate as the
                       %cellular side; satisfied CPEs may still HOST NCRs).
                       %The needy-only subband keeps the forwarded stream
                       %mix small - full-band forwarding under 20+ MU-MIMO
                       %streams was shown to collapse SINR to the mix's
                       %per-stream SIR. 0 = repeater-free FWA phase
params.num_repeater_per_cpe = 2; %NCRs attached per assisted user: TWO independent
                                 %rank-1 amplify-and-forward branches make the
                                 %composite channel rank-2, restoring the UE's
                                 %second spatial stream that a single Rel-18 NCR
                                 %collapses
rep_assist_frac_arr = 0:0.2:1; %swept fraction of each sector's WEAKEST cellular
                                 %UEs (cell edge, by serving-link large-scale gain)
                                 %eligible for repeater assistance; the rest keep the
                                 %pure direct path. 0 disables NCR assistance entirely
                                 %(no attachments), giving the no-repeater baseline
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

%Two static FWA DEMAND CONFIGURATIONS: per-CPE minimum rates by zone
%[apt-block, home, mall, office], anchored to the US FWA plan landscape.
%LOW = the lower end of expected demand per building class; HIGH = the
%upper end. Residential plans span ~50 Mbps entry to ~300 Mbps top
%standard plans (Light Reading,
%https://www.lightreading.com/fixed-wireless-access/verizon-and-t-mobile-fwa-speeds-are-actually-increasing;
%tecknexus, https://tecknexus.com/t-mobile-verizon-and-att-boost-5g-fwa-speeds-and-subscribers/;
%Verizon FWA guide, https://www.verizon.com/home/internet/guides/what-is-fixed-wireless-access-fwa-technology/).
%Business FWA plans run 100/200/400 Mbps (Verizon Business 5G Internet,
%https://www.verizon.com/business/products/internet/5g/): office and
%mall CPEs take 100 at the low end, 200 (the mid business tier) in the
%MID config, and 500 at the high end (top business tier plus margin for
%multi-tenant units). The apartment-block CPE serves an MDU (several
%units behind one CPE): base = the 100 Mbps FCC fixed-broadband
%benchmark tier; mid = 300 (the top standard residential plan as an MDU
%aggregate); high = 500 (assumption, flagged: between the top
%residential plan and the gigabit tier). The MID home rate of 150 Mbps
%is the median delivered US FWA speed (Verizon/AT&T ~150, cited above).
%              zone:   apt   home   mall  office
%Named, swept FWA DEMAND PROFILES: per-zone minimum committed rates
%[apt-block(MDU) home mall office] in Mbps. The profile NAME tags the
%result filenames and the demand_profile CSV column.
fwa_demand_names = {'low','medium','high'};
fwa_demand_configs = 1e6*[100    50    100    100;   %LOW demand profile
                          300   150    200    200;   %MEDIUM demand profile
                          500   300    500    500];  %HIGH demand profile
%% Simulation FR1 setup
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
    for idxnumCPE = 1:length(take_rate_arr) %FWA take-rate sweep over the inventory
    %% CPE locations: numCPE_all per site across THREE zones - apartment
    %top floors, single-family rooftops in the homes belt (FWA's core
    %market), and the strip-mall rooftop at the origin
    %office-park centres: n_officepark per site at opDist, random azimuths
    opC = zeros(M_sites*params.n_officepark,1);
    for idxSite = 1:M_sites
        thOP = 2*pi*rand(params.n_officepark,1);
        opC((idxSite-1)*params.n_officepark+(1:params.n_officepark)) = ...
            siteLocations(idxSite,1) + 1i*siteLocations(idxSite,2) + params.opDist*exp(1i*thOP);
    end
    n_parks = numel(opC);
    %BUILDING INVENTORY per site: apartment blocks (fixed count per
    %complex), single-family homes (residents / household size over the
    %belt area), office buildings, and the shared mall units
    A_home_inv = pi*((params.homeRadius/1000)^2 - (params.aptRadius/1000)^2); %km^2
    inv_apt = params.n_apt_blocks;
    inv_home = round(params.rho_res_home*A_home_inv/params.hh_size);
    inv_op = params.n_officepark*params.n_op_bldg;
    take_rate = take_rate_arr(idxnumCPE);
    n_apt_site = round(take_rate*inv_apt);
    n_home_site = round(take_rate*inv_home);
    %commercial zones use the (near-universal) BUSINESS take rate
    n_op_site = round(params.take_rate_biz*inv_op);
    n_cpe_mall = round(params.take_rate_biz*params.n_mall_units); %shared zone (total)
    n_cpe_op = M_sites*n_op_site;
    numCPE_tot = M_sites*(n_apt_site + n_home_site + n_op_site) + n_cpe_mall;
    numCPE_all = round(numCPE_tot/M_sites); %per-site nominal, for records
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
    %(the demand tier of each CPE follows its zone; set after cpe_zone below)
    params.cpe_zone = [ones(M_sites*n_apt_site,1); 2*ones(M_sites*n_home_site,1); ...
                       3*ones(n_cpe_mall,1); 4*ones(n_cpe_op,1)];
    %per-CPE demand rates are set with the demand configuration below
    for idxActivity = 1:length(activity_arr) %device-activity sweep
        act_scale = activity_arr(idxActivity)/params.activity_ref;
        lambda_road = lambda_UE_road_base*act_scale;
        %% UE construction, per suburban zone. For each zone the GROSS
        %actives split into outdoor pedestrians (pedestrian_ratio, all
        %cellular) and indoor users, of which only (1 - indoor_offload)
        %remain cellular. Every cellular UE carries its zone id, motion
        %region (centre + annulus) and speed, built here explicitly.
        A_apt = pi*((params.aptRadius/1000)^2 - (params.min_dist_2D/1000)^2);  %km^2 per site
        A_home = pi*((params.homeRadius/1000)^2 - (params.aptRadius/1000)^2);  %km^2 per site
        A_mall = pi*(params.mallRadius/1000)^2;                                %km^2, shared
        A_roads = pi*((params.deployRange/1000)^2 - (params.aptRadius/1000)^2);%km^2 per site
        g_apt = ceil(lambda_UE_apt_base*act_scale*A_apt);      %gross actives per site
        g_home = ceil(lambda_UE_home_base*act_scale*A_home);   %gross actives per site
        g_mall = ceil(lambda_UE_mall_base*act_scale*A_mall);   %gross actives, shared zone
        A_park = pi*(params.opRadius/1000)^2;                  %km^2 per office park
        g_park = ceil(lambda_UE_op_base*act_scale*A_park);     %gross actives per park
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
        fadingState = rng; %COMMON RANDOM NUMBERS across the sweep axes below:
        %the state is restored at the top of every pool-size iteration, so
        %every (frac, num_rep, gain, SI) combination replays the identical
        %layer-3 fading draws and sweep comparisons are paired at all three
        %randomness layers (drops still differ across seeds/take rates)
        for idxfrac = 1:length(rep_assist_frac_arr) %repeater-eligibility sweep
        params.rep_assist_frac = rep_assist_frac_arr(idxfrac); %consumed by the rate
            %function's per-sector attachment eligibility and mirrored by the
            %pool selection below; placed after the snapshot generation so
            %every frac value sees the identical drop and motion
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
        for idxnumrep = 1:length(num_rep_arr)  
            rng(fadingState); %replay the same fading draws (see above)
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
                %fixed FWA interference-suppression factor gamma_I
                %(citations at params.SI_cancel_dB above)
                    params.SI_cancel_factor = 10^(0.1*params.SI_cancel_dB);
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
                    %FWA-side NCR assistance: the SAME enabled repeater pool
                    %forwards on the FWA band for the weakest CPEs (same
                    %eligibility fraction, donor-sector lock, and benefit
                    %gate as the cellular side; see the rate function)
                    params.fwa_rep_pool = greedy_order(1:min(params.num_repeater_tot,numCPE_tot));
                    params.rep_P_in_FWA_mW = P_in_rep_mW*params.Band/Band; %constant-PSD
                        %share of the drop-level donor input power that falls in
                        %the FWA band (params.Band = Band_after_cell here)
                    %stationary CPEs: no mobility layer; the fading draw
                    %count matches the cellular sample count (layers 2 x 3)
                    nbrFWADraws = nbrOfSnapshots*nbrOfRealizations;
                    rate_dl = zeros(K,nbrFWADraws);
                    fwaStateFWA = rng; %replayed per demand config below so the
                        %hybrid needy-subband recompute sees the SAME channels
                    for n = 1:nbrFWADraws
                        [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                        %baseline: repeater-free shared full band; the needy
                        %set of each demand config is identified from this
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    save_old_rate = rate_dl;
                    save_old_mean_FWA = mean_rate_dl_FWA;
                    if params.Band > 0
                        %raw FWA packing matrix (skipped if the cell phase
                        %consumed the whole band: all-zero rates carry no
                        %distribution information)
                        writematrix(rate_dl, strcat(packFolder,'/packing_FWA_',num2str(numCPE_all),'CPE_',num2str(params.num_repeater_tot),'rep_',aID,'.csv'));
                    end
                    Band_FWA = params.Band;
                    for idxrmin = 1:size(fwa_demand_configs,1)
                        %per-CPE absolute minimum rates for this config;
                        %the HOME rate labels the config in filenames/CSV
                        params.fwa_demand_rates = fwa_demand_configs(idxrmin, params.cpe_zone).';
                        params.fwa_demand_profile = fwa_demand_names{idxrmin}; %tags filenames/CSV
                        params.r_min_FWA = fwa_demand_configs(idxrmin,2); %internal fallback only
                                              %(computeUtility uses fwa_demand_rates)
                        K_FWA_max = 0;
                        params.Band = Band_FWA;
                        if params.FWA_REPEAT && Band_FWA > 0
                            %HYBRID NCR-assisted scheduling (see the flag):
                            %needy = misses its demand on the shared full band.
                            %The needy KEEP their f_sat share of the shared
                            %band and ADDITIONALLY get the freed (1-f_sat)
                            %slice as a REDUCED MU-MIMO system (channels
                            %sliced to needy CPEs + NCR hosts, so the slice
                            %sees no interference from the satisfied CPEs):
                            %rate_needy = f_sat*old + slice rate
                            rmin_all = params.fwa_demand_rates(:);
                            needy = mean_rate_dl_FWA(1:numCPE_tot) < rmin_all;
                            if any(needy) && any(~needy)
                                %satisfied CPEs keep the smallest shared
                                %subband on which ALL still meet demand
                                f_sat = max(rmin_all(~needy)./mean_rate_dl_FWA(~needy));
                                %reduced node set: needy CPEs + enabled NCR
                                %hosts (hosts only forward on the slice, they
                                %are not scheduled on it)
                                subIdx = union(find(needy), params.fwa_rep_pool(:));
                                needyPos = find(ismember(subIdx, find(needy)));
                                params_sub = params;
                                params_sub.numCPE = numel(subIdx);
                                params_sub.D_FWA = D_FWA(:,subIdx);
                                params_sub.D = params_sub.D_FWA;
                                params_sub.D_cell = [];
                                params_sub.gainOverNoise_lin = params.gainOverNoise_lin(:,subIdx);
                                params_sub.BETA_interUE = params.BETA_interUE(subIdx,subIdx);
                                params_sub.rep_donor_gain = params.rep_donor_gain(:,subIdx);
                                params_sub.set_repeat = find(~ismember(subIdx, find(needy))); %hosts: forward only
                                params_sub.fwa_rep_pool = find(ismember(subIdx, params.fwa_rep_pool));
                                params_sub.rep_assist_frac = params.rep_assist_frac_fwa; %all needy CPEs admitted
                                params_sub.Band = (1-f_sat)*Band_FWA; %the freed slice
                                params_sub.rep_P_in_FWA_mW = P_in_rep_mW(subIdx)*params_sub.Band/Band;
                                rng(fwaStateFWA); %full redraw replays the SAME channels
                                    %as the baseline; the reduced system is a SLICE of them
                                rate_needy = zeros(numel(subIdx),nbrFWADraws);
                                for n = 1:nbrFWADraws
                                    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                                    rate_needy(:,n) = compute_link_rates_MIMO_mmse_wi_repeater(params_sub, channel_dl, channel_est_dl, ...
                                        channel_dl_FWA(:,subIdx,:,:), channel_est_dl_FWA(:,subIdx,:,:), ...
                                        channel_interFWA(subIdx,subIdx,:,:), channel_interFWA_est(subIdx,subIdx,:,:));
                                end
                                rate_dl(needy,:) = f_sat*rate_dl(needy,:) + rate_needy(needyPos,:);
                                rate_dl(~needy,:) = f_sat*rate_dl(~needy,:);
                                mean_rate_dl_FWA = mean(rate_dl,2);
                            end
                        end
                        [cell_util, FWA_util] = computeUtility(params,mean_rate_dl_cell, mean_rate_dl_FWA);
                        K_FWA_init = sum(FWA_util>0);
                        if (params.Band > 0)
                            while params.enable_fwa_cap_realloc && any(mean_rate_dl_FWA > params.r_max_FWA)
                                CPE_idxs = find(mean_rate_dl_FWA > params.r_max_FWA);
                                %scale the cap-exceeders down to (at least)
                                %their OWN tiered demand, freeing the rest
                                rmin_vec = params.fwa_demand_rates;
                                sc = max(rmin_vec(CPE_idxs)./mean_rate_dl_FWA(CPE_idxs));
                                params.Band = params.Band*(1-sc);
                                rate_dl(CPE_idxs,:) = rate_dl(CPE_idxs,:)*sc;
                                mean_rate_dl_FWA(CPE_idxs) = mean(rate_dl(CPE_idxs,:),2);
                                params.set_repeat = [params.set_repeat; CPE_idxs];
                                not_set_repeat = setdiff(1:numCPE_tot,params.set_repeat);
                                for n = 1:nbrFWADraws
                                    [channel_dl, channel_est_dl,channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est] = computePhysicalChannels_sub6_MIMO(params);
                                    if params.FWA_REPEAT
                                        rate_dl(:,n) = compute_link_rates_MIMO_mmse_wi_repeater(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA, channel_interFWA, channel_interFWA_est);
                                    else
                                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl, channel_est_dl, channel_dl_FWA, channel_est_dl_FWA);
                                    end                                              
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
                        result_string = strcat('/results_numFWA_',num2str(numCPE_all),...
                            'CPE_',num2str(activity_arr(idxActivity)),'activity_', ...
                            num2str(deployRange),'deployRange_', ...
                            params.fwa_demand_profile,'demand_', ...
                            num2str(params.num_repeater_tot),'num_repeater_', ...
                            num2str(params.rep_assist_frac),'assistfrac_', ...
                            aID);
                        recording_text_file_string = strcat(rateFolder,result_string,'.csv');
                        fileID = fopen(recording_text_file_string,'w');
                        output_categories = ['numCPE,','activity,','deployRange,','r_min_cell,', ...
                        'demand_profile,','num_rep,','rep_frac,','init_FWA,','max_FWA,','Band,' 'Band_FWA,', 'cell_se,', 'FWA_se\n'];
                        fprintf(fileID,output_categories);
                        formatSpec = '%d,%f,%d,%f,%s,%d,%f,%d,%d,%d,%f,%f,%f\n';
                        fprintf(fileID,formatSpec,numCPE_all,activity_arr(idxActivity),deployRange, ...
                        params.r_min_cell,params.fwa_demand_profile,params.num_repeater_tot,params.rep_assist_frac,K_FWA_init,K_FWA_max,Band,Band_FWA, sum(mean_rate_dl_cell)/(Band - Band_FWA), FWA_se_out);
                        fclose(fileID);
                    end
        end
        end %idxfrac (repeater-eligibility sweep)
    end
    end %idxnumCPE (FWA take-rate sweep)
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)

