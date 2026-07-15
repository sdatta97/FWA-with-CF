%%DEPLOYMENT DIAGRAM GENERATOR: 3D scene (buildings with heights, masts, vegetation) + Blender CSV export
%%Replays the SimulationMain drop construction verbatim (seed 0084 via the
%%default aID path); adjust take_rate_arr/activity_arr overrides below for
%%other operating points. Run from the repo root. Extracted from
%%SimulationMain.m - regenerate the head section if the construction changes.
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
%PER-CLASS CSI FEEDBACK LOOPS: the sounding periodicity is matched to
%each UE class's channel coherence time T_c = 0.423/f_D (Clarke's model;
%Rappaport, "Wireless Communications", Table 5.1) - the heterogeneous-
%coherence training principle (https://arxiv.org/abs/1707.03571). UEs
%whose T_c is shorter than the default 5 ms loop (crossover
%0.423*c/(fc*5ms) = 24.7 km/h at 3.7 GHz) get the FAST loop: NR allows
%SRS/CSI-RS periodicities down to slot level (TS 38.331 SRS-Config /
%TS 38.214 CSI resource periodicities), and adapting them to UE
%Doppler is established practice. In-car UEs (40 km/h, T_c = 3.1 ms):
%1.25 ms loop, staleness rho = J0(2*pi*137*0.625ms) = 0.93; pedestrians
%(3 km/h, T_c = 41 ms): 5 ms loop, rho = 0.993. The residual staleness
%enters as TRUE Jakes aging (estimate formed from the SOUNDED channel,
%truth evolved by rho: Truong & Heath, "Effects of Channel Aging in
%Massive MIMO Systems", https://arxiv.org/abs/1305.6151) and the faster
%sounding is charged to the fast class's prelog (see the prelog block)
params.T_csi_loop_slow = 5e-3;    %default CSI loop (s)
params.T_csi_loop_fast = 1.25e-3; %fast-class CSI loop (s), TS 38.331-configurable
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
%%


%FWA CPE deployment: potential CPE SITES follow the BUILDING INVENTORY
%(one per apartment block, single-family home, office building, and mall
%unit), and the sweep is the TAKE RATE applied to that inventory. US FWA
%serves ~12%% of broadband households in mid-2025, ~15%% by early 2026
%(https://www.telecompetitor.com/fixed-wireless-is-thriving-and-that-could-be-a-problem-report/;
%https://insights.opensignal.com/2025/10/20/the-state-of-us-fwa-what-impact-has-att-internet-airs-launch-had/dt;
%https://www.ericsson.com/en/reports-and-papers/mobility-report/dataforecasts/fwa-outlook),
%so 5-20%% brackets the market. The inventory derives from the resident
%densities and household size below, keeping CPE supply proportional to
%the population it serves.
take_rate_arr = 0.20; %diagram: production take rate FWA take rate applied to the
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
%FWA's market case is strongest where fiber is scarce, i.e. suburbs.
%Two macro sites are deployed CENTRALLY at (+/- ISD/2, 0) along a
%randomly oriented axis, so the deployment region surrounds them.
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
%Two FIXED macro sites (M_sites below)
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
params.rho_res_home = 1200; %residents per km^2, single-family belt
params.activity_ref = 0.05; %REFERENCE busy-hour activity share: the per-zone
                         %gross-active densities below are the values observed
                         %at this share (apartment ring = 5000 residents/km^2 x 5%)
activity_arr = 0.05; %diagram: 5% activity level busy-hour device-activity share; every
                         %zone's gross-active density (and the road density)
                         %scales linearly with activity_arr/activity_ref
lambda_UE_apt_base = 250;  %GROSS active users per km^2 at activity_ref, apartment ring
lambda_UE_home_base = 60;  %GROSS active users per km^2 at activity_ref, single-family belt
lambda_UE_mall_base = 150; %GROSS active users per km^2 at activity_ref, strip-mall parcel
lambda_UE_op_base = 200;   %GROSS active users per km^2 at activity_ref, office-park
                         %parcels (sparse evening staff/visitors; the parks' main
                         %role in the busy hour is hosting CPEs/NCRs)
params.incar_share = 0.20; %fraction of the gross active population that is
                         %outdoor IN CARS: the 3GPP evaluation user-
                         %distribution convention "80% indoor, 20% of users
                         %outdoor in cars" (TR 38.913 urban coverage
                         %attributes; same split in the TR 38.901 clause 7.8
                         %calibration). Replaces the legacy independent road
                         %density, which at ISD 1732 made the vehicular class
                         %half the cellular population
params.pedestrian_ratio = 0.2; %fraction of each zone's GROSS actives outdoors on foot
params.indoor_offload = 0.8;   %fraction of indoor actives on home/venue broadband
                               %instead of cellular (Cisco VNI offload + US home
                               %internet penetration, see above)
params.loss_pc_cell = 5/100; %enforce r_min_cell at the 5th-percentile cell UE
       %(ITU-R M.2410 / TR 38.913 5th-percentile user-SE KPI convention) -
       %the band-setting UE then sits at the noise-limited edge where NCR
       %assistance can act. The FCC BDC mobile-coverage rule would support
       %10% (>= 90% cell-edge coverage probability at 35/3 Mbps, >= 50%
       %loading: https://help.bdc.fcc.gov/hc/en-us/articles/6047425308187),
       %but at the 10th percentile the setter is interference-limited and
       %decoupled from the repeater mechanism; recorded in the loss_pc
       %naming-registry entry either way
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
params.HW_IMPAIRMENTS = 0;  % 1 = hardware impairments on, 0 = ideal hardware
                            %(recorded via the HWI naming-registry entry)
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
params.fwa_pilot_contam = 0; %DORMANT (=1 to study): PILOT CONTAMINATION for
                            %the MU-MIMO FWA CPEs, whose downlink precoding CSI
                            %comes from TDD-reciprocal UPLINK SRS sounding
                            %(TS 38.214 clause 6.2.1); cellular UEs keep the
                            %per-link CSI-RS estimation below
params.tau_p_fwa = 48;      %orthogonal SRS resources per site for CPE sounding:
                            %comb-4 (K_TC = 4 comb offsets) x 12 cyclic shifts
                            %in one SRS symbol of the 5 ms sounding period
                            %(TS 38.211 clause 6.4.1.4); CPEs beyond tau_p per
                            %site REUSE pilots -> full intra-site collision
params.srs_root_leak_dB = -29; %cross-site SRS leakage power: neighbouring sites
                            %use different Zadoff-Chu sequence groups (TS 38.211
                            %clause 6.4.1.4.1 group planning), whose
                            %cross-correlation power ~ 1/M_sc with
                            %M_sc = 273 PRB x 12 / 4 = 819 comb subcarriers at
                            %100 MHz -> -29 dB; cyclic shifts do not separate
                            %ACROSS roots, so all co-comb-offset resources of
                            %the other site leak at this level
params.csi_rs_offset_dB = 0;% CSI-RS EPRE relative to PDSCH EPRE (powerControlOffset,
                            % TS 38.214; standardised range -8..15 dB, default 0)
%Network-controlled repeater (NCR) amplification, per 3GPP: fixed
%amplification gain of 90 dB subject to the maximum output power
%constraint of 33 dBm (evaluation assumptions of TR 38.867, "Study on NR
%network-controlled repeaters"; RF repeater requirements in TS 38.106).
%The cap is applied per repeater at its actual operating point below, so
%the effective gain is min(G_max, P_max/P_in)
%PILOT-SIGNALLING PRELOG OVERHEADS, applied only under IMPERFECT_CSI
%(perfect CSI is genie-aided: preLogFactor = 1, no pilots). 30 kHz SCS:
%the 5 ms sounding period spans 10 slots x 14 = 140 OFDM symbols.
% - cellular OFDMA UEs (PER UE, set at drop level with the per-class
%   CSI loops): PDSCH DMRS type 1, single symbol per slot (TS 38.211
%   clause 7.4.1.1) + 64-port CSI-RS as 2x 32-port resources, density 1
%   RE/port/RB (TS 38.211 clause 7.4.1.5) per that UE's sounding period
%   - the 1.25 ms fast class pays 4x the CSI-RS overhead;
% - TDD MU-MIMO FWA CPEs: double-symbol DMRS (extended orthogonal port
%   set for MU-MIMO, TS 38.211 Table 7.4.1.1.2-2) + the SRS sounding
%   symbols, ceil(CPEs-per-site/tau_p_fwa) of the 140 symbols
%   (set per drop after the SRS pilot plan; see preLogFactor_fwa there).
params.preLogFactor = 1; %per-phase value: the cellular phase assigns the
                         %per-UE preLogFactor_cell vector (drop level), the
                         %FWA phase the scalar preLogFactor_fwa
params.repeat_gain = 90; %FIXED NCR amplification gain G_max (dB), TR 38.867 (not swept)
params.rep_max_pow_dBm = 33; %NCR maximum output power (dBm), TR 38.867 / TS 38.106
params.ncr_zf_donor = 1; %1 = interference-aware (max-SINR / regularized
                         %zero-forcing) donor receive combining at the NCR,
                         %suppressing the co-channel sites at the donor input
                         %- the SIR ceiling of the repeated path; 0 = matched
                         %combining only (ablation)
params.ncr_rank = 1; %NCR forwarded spatial layers: 1 = Rel-18 single-chain
                     %AF (Carvalho et al., IEEE Comm. Standards Mag. 2026,
                     %https://arxiv.org/pdf/2403.09601); params.N_UE_cell =
                     %beyond-Rel-18 MIMO-repeater what-if (details in
                     %compute_link_rates_OFDM_wi_repeater.m)
num_rep_arr = 0:1:12; %total (0 = no-NCR baseline: sources the plain-cellular
                       %packing matrices and the pool-sweep reference) enabled repeaters; greedy max-coverage attachment
params.ncr_benefit_gate = 1; %per-UE NCR ON/OFF side control (TR 38.867); =0 only for ablation
                     %is per donor sector, saturating near one per sector (6)
params.FWA_REPEAT = 0; %DISABLED for now (avoids confounding the cellular-side
                       %NCR study). 1 = HYBRID FWA scheduling with NCR assistance: CPEs
                       %missing their demand on the shared band ("needy") get
                       %the freed remainder of the satisfied CPEs' subband
                       %(shrunk so every satisfied CPE still meets demand) as
                       %a dedicated NCR-aided slice ON TOP of their f_sat
                       %share; scheme and rationale documented in
                       %compute_link_rates_MIMO_mmse_wi_repeater.m and at the
                       %recompute in the demand loop. 0 = repeater-free
params.num_repeater_per_cpe = 1; %NCRs attached per assisted user (single
                                 %rank-1 amplify-and-forward branch; 2 restores
                                 %the UE's second spatial stream via two
                                 %independent branches)
%Three layers of randomness, averaged per SLURM array task:
%  (1) the seed aID fixes the DROP: gNB/CPE/UE locations and the initial
%      large-scale state;
%  (2) nbrOfSnapshots units of MOBILITY: the cellular UEs advance
%      dt_snap along their trajectories and the large-scale state
%      (shadowing, LOS, association) evolves consistently. Motion model:
%      straight-line constant speed with a random heading per user
%      (Random Direction model: Camp, Boleng & Davies, WCMC 2002,
%      https://onlinelibrary.wiley.com/doi/10.1002/wcm.72), the 3GPP /
%      ITU-R M.2412 evaluation convention; in-car road UEs at 40 km/h in the
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

%Named, swept FWA DEMAND PROFILES: per-zone minimum committed rates
%[apt-block(MDU) home mall office] in Mbps; the profile NAME tags the
%result filenames and the demand_profile CSV column. Anchors: US
%residential plans span ~50 Mbps entry to ~300 Mbps top standard tiers,
%~150 Mbps median delivered (Light Reading,
%https://www.lightreading.com/fixed-wireless-access/verizon-and-t-mobile-fwa-speeds-are-actually-increasing;
%tecknexus, https://tecknexus.com/t-mobile-verizon-and-att-boost-5g-fwa-speeds-and-subscribers/;
%Verizon FWA guide, https://www.verizon.com/home/internet/guides/what-is-fixed-wireless-access-fwa-technology/).
%Business plans run 100/200/400 Mbps (Verizon Business 5G Internet,
%https://www.verizon.com/business/products/internet/5g/): offices/malls
%take 100/200/500 across the profiles (HIGH = top tier + multi-tenant
%margin). The apartment-block CPE serves an MDU: 100 (FCC fixed-broadband
%benchmark) / 300 (top residential plan as aggregate) / 500 (assumption,
%flagged).
fwa_demand_names = {'low','medium','high'};
fwa_demand_configs = 1e6*[100    50    100    100;   %LOW demand profile
                          300   150    200    200;   %MEDIUM demand profile
                          500   300    500    500];  %HIGH demand profile
%% AUTOMATED OUTPUT NAMING (see sweepNaming.m): every sweepable
%parameter is registered here; axes with numel > 1 are detected as
%swept and become CSV table columns, scalars become constants whose
%NAMES form the result folder and whose VALUES form the file name.
%Register any new sweep axis in BOTH lists (and in currVals at the
%writer) - the check_naming.m scratchpad script verifies consistency
%without running the simulation.
sweepRegNames = {'take_rate','activity','num_rep','demand_profile', ...
    'rep_gain','SI_cancel','r_min_cell_Mbps','loss_pc','ISD_m','Band_MHz','HWI', ...
    'FWArep','rank','Krep','gate','contam'};
sweepRegVals = {take_rate_arr, activity_arr, num_rep_arr, ...
    fwa_demand_names, params.repeat_gain, params.SI_cancel_dB, ...
    params.r_min_cell/1e6, params.loss_pc_cell, params.ISD, Band/1e6, params.HW_IMPAIRMENTS, ...
    params.FWA_REPEAT, params.ncr_rank, ...
    params.num_repeater_per_cpe, params.ncr_benefit_gate, params.fwa_pilot_contam};
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
        %% UE construction, per suburban zone. For each zone the GROSS
        %actives split into outdoor pedestrians (pedestrian_ratio, all
        %cellular) and indoor users, of which only (1 - indoor_offload)
        %remain cellular. Every cellular UE carries its zone id, motion
        %region (centre + annulus) and speed, built here explicitly.
        A_apt = pi*((params.aptRadius/1000)^2 - (params.min_dist_2D/1000)^2);  %km^2 per site
        A_home = pi*((params.homeRadius/1000)^2 - (params.aptRadius/1000)^2);  %km^2 per site
        A_mall = pi*(params.mallRadius/1000)^2;                                %km^2, shared
        g_apt = ceil(lambda_UE_apt_base*act_scale*A_apt);      %gross actives per site
        g_home = ceil(lambda_UE_home_base*act_scale*A_home);   %gross actives per site
        g_mall = ceil(lambda_UE_mall_base*act_scale*A_mall);   %gross actives, shared zone
        A_park = pi*(params.opRadius/1000)^2;                  %km^2 per office park
        g_park = ceil(lambda_UE_op_base*act_scale*A_park);     %gross actives per park
        %in-car actives per site: cars are incar_share of the TOTAL gross
        %active population (TR 38.913 convention, see the param), i.e.
        %share/(1-share) of the non-car zone actives; they scale with the
        %activity sweep through the zone gross counts
        n_car = ceil(params.incar_share/(1-params.incar_share)* ...
            (g_apt + g_home + g_mall/M_sites + params.n_officepark*g_park)); %per site (all cellular)
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


    %% --- 3D deployment diagram: buildings, masts, vegetation ---
    rng(7); %visual-only randomness (tree placement, unmodeled heights)
    fig = figure('Position',[80 60 1000 760], 'Color','w');
    hold on
    %ground disc per site + vegetation annulus
    th = linspace(0,2*pi,120);
    for s = 1:M_sites
        cx = siteLocations(s,1); cy = siteLocations(s,2);
        patch(cx+params.deployRange*cos(th), cy+params.deployRange*sin(th), ...
              zeros(size(th)), [0.93 0.92 0.88], 'EdgeColor','none');
    end
    for s = 1:M_sites
        cx = siteLocations(s,1); cy = siteLocations(s,2);
        rv = [params.homeRadius params.deployRange];
        xx = [cx+rv(1)*cos(th), cx+rv(2)*cos(fliplr(th))];
        yy = [cy+rv(1)*sin(th), cy+rv(2)*sin(fliplr(th))];
        patch(xx, yy, 0.5*ones(size(xx)), [0.72 0.85 0.66], 'EdgeColor','none','FaceAlpha',0.7);
    end
    hVeg = patch(NaN, NaN, [0.72 0.85 0.66], 'EdgeColor','none'); %legend proxies
    hApt = patch(NaN, NaN, [0.35 0.55 0.80], 'EdgeColor',[.35 .35 .35]);
    hHom = patch(NaN, NaN, [0.85 0.80 0.62], 'EdgeColor',[.35 .35 .35]);
    hMal = patch(NaN, NaN, [0.66 0.44 0.66], 'EdgeColor',[.35 .35 .35]);
    hOff = patch(NaN, NaN, [0.90 0.55 0.30], 'EdgeColor',[.35 .35 .35]);
    hCpe = patch(NaN, NaN, [0.05 0.15 0.5], 'EdgeColor','none');
    hTre = patch(NaN, NaN, [0.23 0.49 0.25], 'EdgeColor','none');
    hCar = patch(NaN, NaN, [0.75 0.10 0.15], 'EdgeColor','none');
    %helper: extruded box
    drawBox = @(x,y,w,h,c) patch( ...
        'Vertices',[x-w/2 y-w/2 0; x+w/2 y-w/2 0; x+w/2 y+w/2 0; x-w/2 y+w/2 0; ...
                    x-w/2 y-w/2 h; x+w/2 y-w/2 h; x+w/2 y+w/2 h; x-w/2 y+w/2 h], ...
        'Faces',[1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 5 6 7 8], ...
        'FaceColor',c,'EdgeColor',[.35 .35 .35],'LineWidth',0.2,'FaceAlpha',1);
    %buildings at CPE sites, heights per zone floor conventions (3 m/floor)
    zc = params.cpe_zone;
    for k2 = 1:size(CPE_locations,1)
        x = CPE_locations(k2,1); y = CPE_locations(k2,2);
        switch zc(k2)
            case 1, fl = randi([2 4]);  w = 28; col = [0.35 0.55 0.80]; %apartment block
            case 2, fl = randi([1 2]);  w = 12; col = [0.85 0.80 0.62]; %single-family home
            case 3, fl = 1;             w = 30; col = [0.66 0.44 0.66]; %mall unit
            case 4, fl = randi([2 3]);  w = 24; col = [0.90 0.55 0.30]; %office building
        end
        drawBox(x, y, w, 3*fl, col);
        %rooftop CPE: small dark cap on the roof instead of a floating marker
        patch('Vertices',[x-3 y-3 3*fl+0.5; x+3 y-3 3*fl+0.5; x+3 y+3 3*fl+0.5; x-3 y+3 3*fl+0.5], ...
              'Faces',[1 2 3 4],'FaceColor',[0.05 0.15 0.5],'EdgeColor','none');
    end
    %vegetation: tree glyphs (cone canopies) in the buffer annulus
    for s = 1:M_sites
        cx = siteLocations(s,1); cy = siteLocations(s,2);
        nT = 130;
        rr = sqrt(params.homeRadius^2 + (params.deployRange^2-params.homeRadius^2)*rand(nT,1));
        aa = 2*pi*rand(nT,1);
        for k2 = 1:nT
            tx = cx+rr(k2)*cos(aa(k2)); ty = cy+rr(k2)*sin(aa(k2));
            hT = 8 + 6*rand;
            [cxx,cyy,czz] = cylinder([6 0],8);
            surf(tx+cxx, ty+cyy, hT*0.45 + czz*hT*0.75, 'FaceColor',[0.23 0.49 0.25],'EdgeColor','none');
        end
    end
    %gNB masts (35 m) with panel heads
    for s = 1:M_sites
        cx = siteLocations(s,1); cy = siteLocations(s,2);
        [mx,my,mz] = cylinder(4,10);
        surf(cx+mx, cy+my, mz*params.ht_bs, 'FaceColor',[0.25 0.25 0.28],'EdgeColor','none');
        drawBox(cx, cy, 22, 4, [0.15 0.15 0.18]); %equipment shelter footprint
        hGnb = plot3(cx, cy, params.ht_bs+3, '^', 'MarkerSize',6, 'Color','k', 'MarkerFaceColor','k');
    end
    %active UEs: pedestrians/indoor (1.5 m markers), cars (small boxes)
    fast = v_move_ue > 5; ped = ~fast & ~isIndoorUE; ind = ~fast & isIndoorUE;
    hPed = plot3(UE_locations(ped,1), UE_locations(ped,2), 4*ones(sum(ped),1), 'o', 'Color',[0.15 0.4 0.08],'MarkerFaceColor',[0.35 0.68 0.2],'MarkerSize',3);
    hInd = plot3(UE_locations(ind,1), UE_locations(ind,2), 4*ones(sum(ind),1), 'o', 'Color',[0.15 0.4 0.08],'MarkerSize',3);
    carIdx = find(fast);
    for k2 = carIdx'
        drawBox(UE_locations(k2,1), UE_locations(k2,2), 10, 3.5, [0.75 0.10 0.15]);
    end
    %view, lighting, vertical exaggeration (heights x ~6 for visibility)
    lg = legend([hGnb hApt hHom hMal hOff hCpe hPed hInd hCar hTre hVeg], ...
        {'gNB site, 35 m mast (3 sectors)','Apartment block (2-4 floors)', ...
         'Single-family home (1-2 floors)','Strip-mall unit','Office building (2-3 floors)', ...
         'Rooftop FWA CPE','Active pedestrian UE (3 km/h)','Active indoor UE (3 km/h)', ...
         'Active in-car UE (40 km/h)','Vegetation (dense clutter)','Vegetation buffer zone'}, ...
        'Location','northwest','FontSize',8);
    lg.Box = 'off';
    axis off; grid off
    daspect([1 1 0.15]); view(-38, 32);
    camlight('headlight'); lighting gouraud; material dull
    exportgraphics(fig, 'plots/deployment_3d_5pct.png', 'Resolution', 240);
    savefig(fig, 'plots/deployment_3d_5pct.fig');
    fprintf('3D diagram saved (%d buildings, %d active UEs)\n', size(CPE_locations,1), size(UE_locations,1));
    %also export the scene geometry for the Blender pipeline
    fid = fopen('plots/deployment_scene.csv','w');
    fprintf(fid,'type,x,y,height,width\n');
    for s = 1:M_sites, fprintf(fid,'gnb,%.1f,%.1f,%.1f,8\n', siteLocations(s,1), siteLocations(s,2), params.ht_bs); end
    for k2 = 1:size(CPE_locations,1)
        tt = {'apt','home','mall','office'};
        fprintf(fid,'%s,%.1f,%.1f,0,0\n', tt{zc(k2)}, CPE_locations(k2,1), CPE_locations(k2,2));
    end
    for k2 = 1:size(UE_locations,1)
        if fast(k2), uet='car'; elseif isIndoorUE(k2), uet='indoor'; else, uet='ped'; end
        fprintf(fid,'%s,%.1f,%.1f,1.5,0\n', uet, UE_locations(k2,1), UE_locations(k2,2));
    end
    fprintf(fid,'veg,%.1f,%.1f,%.1f,%.1f\n', 0, 0, params.homeRadius, params.deployRange);
    fclose(fid);
    end %idxActivity
    end %idxnumCPE
