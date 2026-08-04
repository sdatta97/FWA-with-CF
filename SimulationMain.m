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
% - mobile 5G-NR coverage benchmark 35/3 Mbps DL/UL (the r_min_cell anchor)
params.r_min_cell = 35e6; %cellular DL floor at the 5th percentile: the FCC
                          %mobile 5G-NR coverage benchmark of 35/3 Mbps DL/UL
                          %(FCC 24-27); attainable here because the FWA/Wi-Fi
                          %offload carries the indoor demand
%%
%%


%FWA CPE deployment: potential CPE SITES follow the BUILDING INVENTORY
%(one per townhouse, single-family home, office building, and mall
%unit), and the sweep is the TAKE RATE applied to that inventory. US FWA
%serves ~12%% of broadband households in mid-2025, ~15%% by early 2026
%(https://www.telecompetitor.com/fixed-wireless-is-thriving-and-that-could-be-a-problem-report/;
%https://insights.opensignal.com/2025/10/20/the-state-of-us-fwa-what-impact-has-att-internet-airs-launch-had/dt;
%https://www.ericsson.com/en/reports-and-papers/mobility-report/dataforecasts/fwa-outlook),
%and Ericsson forecasts ~18%% of US fixed broadband by 2030, so 5-20%%
%brackets the market. The inventory derives from the resident densities
%and household size below, keeping CPE supply proportional to the
%population it serves.
take_rate_arr = 0.20; %DEFAULT residential FWA take rate applied to the
                      %townhouse and single-family inventory: the top of
                      %the 5-20% market bracket (Ericsson 2030 forecast
                      %~18% of fixed broadband; iGR addressable market
                      %25.4% of households). Sweep 0.05:0.05:0.20 for
                      %sensitivity studies.
params.take_rate_biz = 0.95; %BUSINESS take rate for the commercial inventory
       %(office buildings, mall units): "nearly all businesses (95%)
       %report having at least one broadband Internet connection" (FCC
       %Business Broadband Capability Survey,
       %https://www.fcc.gov/document/business-broadband-capability-survey-results;
       %see also https://advocacy.sba.gov/2026/01/13/issue-brief-no-23-small-businesses-access-to-broadband-internet/);
       %in this FWA-served deployment those connections are the CPEs
params.hh_size = 2.5;  %persons per household (US Census ~2.5)
params.rho_res_ring = 2500; %residents per km^2, inner TOWNHOUSE ring
       %(one potential CPE per townhouse, individual service - no MDU
       %aggregation). US suburbs hold few apartment blocks: ~62% of US
       %housing is single-family detached and many suburbs exceed 90%
       %(Census/EIA, https://www.eia.gov/consumption/residential/data/2024/hc/pdf/HC9.1_2024.pdf;
       %https://constructioncoverage.com/research/cities-with-the-most-single-family-homes).
       %Townhouse tracts run 10-20 units/acre NET (ULI/MRSC,
       %https://urbanland.uli.org/development-business/thin-micro-townhouses-optimize-density;
       %https://mrsc.org/stay-informed/mrsc-insight/april-2017/visualizing-compatible-density);
       %the GROSS ring average here (streets, yards, commons) is 1,000
       %units/km^2 = ~4 u/ac -> 2,500 residents/km^2, between the inner-
       %suburb ~1,550/km^2 mean and net townhouse-tract density
       %(https://www.jchs.harvard.edu/blog/the-implications-of-different-suburban-definitions)
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
% - townhouse ring [min_dist_2D, aptRadius] around each site: the denser
%   inner neighbourhood of ATTACHED single-family homes (townhouses,
%   2-3 floors, FWA CPEs on top floors) - suburban housing is
%   differently sized single-family stock, not apartment blocks;
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
params.aptRadius = 175;  %townhouse-ring outer radius (m); the apt_* field
                         %names are kept for schema stability, zone 1 IS
                         %the townhouse ring throughout
params.homeRadius = 450; %single-family belt outer radius (m)
params.mallRadius = 75;  %strip-mall disc radius (m), centred at the origin
%OFFICE PARKS: smaller commercial clusters on the arterial roads in the
%outer belt - low-rise office buildings with FWA CPEs on their rooftops,
%whose co-located NCRs sit close to the noise-limited road users at the
%coverage-limited cell edge (the CPE-hosted deployment logic of the
%paper: repeaters go where subscriber equipment already is)
params.n_officepark = 1; %office parks per site (one modest park is the
                         %defensible suburban norm; was 2)
params.opRadius = 60;    %office-park parcel radius (m)
params.opDist = 600;     %office-park centre distance from its site (m),
                         %on the roads between the homes belt and the edge
params.op_floors = [2 3]; %low-rise office buildings, floors
params.apt_floors = [2 3];  %townhouses, floors (zone 1; field name kept)
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
params.N_UE_cell = 2;
params.hr = 1.5; %outdoor UT height (m)
%CPE mounting height is drawn per CPE in generateSetup.m: top floor of
%its townhouse (params.apt_floors) or the strip-mall rooftop
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
%on ANY access network, per suburban zone. The inner townhouse ring
%holds ~2,500 residents/km^2 (see rho_res_ring: gross ~4 units/acre
%against the 10-20 u/ac NET of townhouse tracts) and single-family
%suburbs ~1,000-1,500 residents/km^2 (Demographia,
%http://demographia.com/db-intlsub.htm; median US suburban tracts
%~700-770/km^2 with inner suburbs ~1,550/km^2,
%https://www.jchs.harvard.edu/blog/the-implications-of-different-suburban-definitions),
%a 2:1 gradient that spreads population evenly instead of concentrating
%it in apartment blocks. A ~5% busy-hour active share - the top of the
%typical 3-5% suburban range, reflecting high-data-rate evening
%applications - gives ~125 and ~60 gross active users/km^2. The strip
%mall carries a modest evening crowd (~150/km^2 over its small parcel).
%Suburban roads carry ~4-15 active in-car users/km^2. The resulting
%~10 cellular-active UEs per sector matches the 10-per-TRxP full-buffer
%evaluation convention of ITU-R M.2412 / 3GPP TR 38.913; at the DEFAULT
%10% activity the ~20 per TRxP stays within the M.2412 framework, whose
%evaluation-configuration Note 6 states higher user densities (e.g. 20
%per TRxP) are not precluded.
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
% - OUTDOOR PEDESTRIANS (3 km/h, no O2I, 1.5 m) in the townhouse ring,
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
                         %at this share (townhouse ring = 2500 residents/km^2 x 5%)
activity_arr = 0.05; %DEFAULT busy-hour device-activity share: the
                         %ITU-R-anchored reference (~10 cellular actives per
                         %TRxP, M.2412 evaluation convention; top of the
                         %typical 3-5% suburban range). Every zone's
                         %gross-active density and the road density scale
                         %linearly with activity_arr/activity_ref; sweep
                         %0.025:0.025:0.1 for sensitivity studies
lambda_UE_apt_base = 125;  %GROSS active users per km^2 at activity_ref,
                           %townhouse ring (rho_res_ring x activity_ref)
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
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10386719/). Set to the
%Tarana ngFWA figure (45 dB): gamma_I applies ONLY to OTHER-SITE
%interferers (see compute_link_rates_MIMO_mmse.m) - same-sector streams
%are suppressed by the L-MMSE precoder itself and co-site sectors arrive
%from the serving mast's direction, so the deep null is spent entirely
%on the one interferer class the CPE can actually steer against
params.SI_cancel_dB = -45;
params.gI_site_aware = 1; %1 = gamma_I nulls OTHER-SITE interferers only
                          %(L-MMSE handles same-sector streams; co-site
                          %sectors arrive from the serving direction);
                          %0 = legacy blanket suppression (ablation).
                          %REGISTERED as a naming constant so this
                          %model generation never shares a result folder
                          %with pre-site-aware campaigns
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
%[townhouse home mall office] in Mbps; the profile NAME tags the
%result filenames and the demand_profile CSV column. Anchors: US
%residential plans span ~50 Mbps entry to ~300 Mbps top standard tiers,
%~150 Mbps median delivered (Light Reading,
%https://www.lightreading.com/fixed-wireless-access/verizon-and-t-mobile-fwa-speeds-are-actually-increasing;
%tecknexus, https://tecknexus.com/t-mobile-verizon-and-att-boost-5g-fwa-speeds-and-subscribers/;
%Verizon FWA guide, https://www.verizon.com/home/internet/guides/what-is-fixed-wireless-access-fwa-technology/).
%Business plans run 100/200/400 Mbps (Verizon Business 5G Internet,
%https://www.verizon.com/business/products/internet/5g/): offices/malls
%take 100/200/500 across the profiles (HIGH = top tier + multi-tenant
%margin). Townhouse CPEs serve ONE household each (no MDU aggregation),
%so both residential zones share the per-home tiers.
fwa_demand_names = {'low'}; %DEFAULT: the LOW profile - the band-limited
%tier where subscriber growth tracks the NCR-freed bandwidth ~1:1
%(medium/high are SINR-limited and sublinear). For tier sensitivity
%restore {'low','medium','high'} with rows [150 150 200 200;
%300 300 500 500] below the low row.
fwa_demand_configs = 1e6*[ 50    50    100    100];  %LOW demand profile
%% AUTOMATED OUTPUT NAMING (see sweepNaming.m): every sweepable
%parameter is registered here; axes with numel > 1 are detected as
%swept and become CSV table columns, scalars become constants whose
%NAMES form the result folder and whose VALUES form the file name.
%Register any new sweep axis in BOTH lists (and in currVals at the
%writer) - the check_naming.m scratchpad script verifies consistency
%without running the simulation.
%'logdet' is a pure GENERATION MARKER, not a model switch: rates are the
%MIMO log-det capacity since 8d957d9 (per-stream SINR before). It keys
%the result folders (and the packing/SE folder suffixes below) so
%log-det campaign files can never pool with per-stream-era data - the
%file tags do not encode the rate model.
logdet_generation = 1;
sweepRegNames = {'take_rate','activity','num_rep','demand_profile', ...
    'rep_gain','SI_cancel','r_min_cell_Mbps','loss_pc','ISD_m','Band_MHz','HWI', ...
    'rank','Krep','gate','contam','gIsite','logdet'};
sweepRegVals = {take_rate_arr, activity_arr, num_rep_arr, ...
    fwa_demand_names, params.repeat_gain, params.SI_cancel_dB, ...
    params.r_min_cell/1e6, params.loss_pc_cell, params.ISD, Band/1e6, params.HW_IMPAIRMENTS, ...
    params.ncr_rank, ...
    params.num_repeater_per_cpe, params.ncr_benefit_gate, params.fwa_pilot_contam, ...
    params.gI_site_aware, logdet_generation};
naming = sweepNaming(sweepRegNames, sweepRegVals, aID);
%a rerun of the same seed starts its table fresh (rows accumulate across
%the sweep loops within one run, never across runs)
if isfile(fullfile('resultData', naming.folder, naming.file))
    delete(fullfile('resultData', naming.folder, naming.file));
end
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
    %% CPE locations: numCPE_all per site across the zones - townhouse
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
    %BUILDING INVENTORY per site: townhouses and single-family homes
    %(residents / household size over each zone's area), office
    %buildings, and the shared mall units
    A_ring_inv = pi*((params.aptRadius/1000)^2 - (params.min_dist_2D/1000)^2); %km^2
    A_home_inv = pi*((params.homeRadius/1000)^2 - (params.aptRadius/1000)^2); %km^2
    inv_apt = round(params.rho_res_ring*A_ring_inv/params.hh_size);
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
    %CPE zone ids: 1 townhouse, 2 single-family home, 3 strip mall, 4 office park
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
        [gainOverNoisedB,gainOverNoisedB_ue,R_gNB,R_cpe,R_interue,R_ue,D_FWA,D_cell,repDonorGaindB,mobState] = generateSetup(params,str2double(aID));
        params.rep_donor_gain = db2pow(repDonorGaindB); %donor-side repeater sector antenna gain (linear power, M_sectors x K_FWA)
        params.gainOverNoise_lin = db2pow(gainOverNoisedB); %BS-to-user large-scale gain (linear), single source of truth
        %SRS pilot plan for the FWA CPEs (drop-level): round-robin over the
        %tau_p orthogonal resources within each site, DETERMINISTIC in the
        %CPE index (consumes no randomness, so the drop structure and the
        %CRN pairing are untouched); the two sites use different sequence
        %roots and share resource numbering (see the params above)
        servSector_cpe = zeros(numCPE_tot,1);
        for q = 1:numCPE_tot
            servSector_cpe(q) = find(D_FWA(:,q)==1,1);
        end
        params.site_of_cpe = ceil(servSector_cpe/S);
        pilot_fwa = zeros(numCPE_tot,1);
        for s = 1:M_sites
            members = find(params.site_of_cpe == s);
            pilot_fwa(members) = mod((0:numel(members)-1)', params.tau_p_fwa) + 1;
        end
        params.pilot_fwa = pilot_fwa;
        %FWA prelog: SRS symbols scale with the per-site CPE count (the
        %pilot pool is tau_p_fwa resources per SRS symbol; see above)
        if params.IMPERFECT_CSI
            n_srs_symb = ceil(max(accumarray(params.site_of_cpe,1))/params.tau_p_fwa);
            params.preLogFactor_fwa = 1 - 2/14 - n_srs_symb/140;
        else
            params.preLogFactor_fwa = 1;
        end
        %Total downlink power arriving at each repeater's donor panel (mW),
        %summed over all sector BSs: sets the operating point of the NCR
        %output-power cap. gainOverNoise_lin is normalized by the full-band
        %noise power, which is re-applied here
        noise_mW = 10^(0.1*(-174 + 10*log10(Band) + params.noiseFigure));
        P_in_rep_mW = params.rho_tot*noise_mW*sum(params.gainOverNoise_lin(:,1:numCPE_tot).*params.rep_donor_gain,1)';
        %Per-UE channel aging factor from the per-UE speeds built above
        %(pedestrians and indoor users 3 km/h, in-car 40 km/h)
        %per-UE CSI loop selection and TRUE staleness factor (see the
        %per-class CSI block in the params section for the physics and
        %references): fast loop when T_c < the default loop; the average
        %sounding-to-use lag is half the loop period
        f_D_ue = v_move_ue*params.fc/params.c;
        fastCSI_ue = 0.423./max(f_D_ue,eps) < params.T_csi_loop_slow;
        T_loop_ue = params.T_csi_loop_slow*ones(size(v_move_ue));
        T_loop_ue(fastCSI_ue) = params.T_csi_loop_fast;
        params.mob_rho = besselj(0,2*pi*(T_loop_ue/2).*f_D_ue);
        %per-UE cellular prelog: the fast class pays proportionally more
        %CSI-RS overhead (64 ports, density 1 RE/port/RB, refreshed every
        %T_loop_ue; UE-dedicated beamformed NZP-CSI-RS resources)
        if params.IMPERFECT_CSI
            symb_per_loop = 140*(T_loop_ue/params.T_csi_loop_slow);
            params.preLogFactor_cell = 1 - 1/14 - 64./(12*symb_per_loop);
        else
            params.preLogFactor_cell = ones(size(v_move_ue));
        end
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
                [snapGain{t_snap},snapGainUE{t_snap},snapR_gNB{t_snap},~,snapR_int{t_snap},~,~,snapDcell{t_snap},~,mobState] = generateSetup(params,0,mobState);
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
        %% Enabled-repeater pool order: GREEDY MAXIMUM COVERAGE over the
        %cellular UEs. EVERY served UE is eligible for assistance (the
        %per-UE benefit gate in the rate function performs the actual
        %selection); a repeater can cover a UE iff its donor sector is
        %the UE's serving sector and the donor x service link metric is
        %positive.
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
        elig_idx = find(sum(D_cell,1) > 0); %every served UE
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
                params.preLogFactor = params.preLogFactor_cell; %cellular-phase prelog
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
                %packing/SE folders carry the _logdet generation suffix
                %like the results folders: the file tags (CPE/act/rep/seed)
                %do not encode the rate/interference model, so per-stream-
                %era matrices must never share a folder with these
                packFolder = 'resultData/FWA_packing_analysis_logdet';
                if not(isfolder(packFolder))
                    mkdir(packFolder)
                end
                if params.CELL_REPEAT
                    cellTag = 'packing_cell_rep_';
                else
                    cellTag = 'packing_cell_plain_';
                end
                if num_rep_arr(idxnumrep) == 0 %packing matrices ONLY at the
                    %no-NCR baseline (the packing figure's configuration; the
                    %activity token keeps the sweep's operating points apart)
                    writematrix(rate_dl, strcat(packFolder,'/',cellTag,num2str(numCPE_all),'CPE_',num2str(activity_arr(idxActivity)),'act_',num2str(params.num_repeater_tot),'rep_',aID,'.csv'));
                end
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
                    params.preLogFactor = params.preLogFactor_fwa; %FWA-phase prelog
                    params.numUE = 0;
                    params.numCPE = numCPE_tot;
                    params.CPE_locations = CPE_locations;
                    params.D_FWA = D_FWA;
                    params.R_gNB = R_gNB(:,:,:,1:numCPE_tot);
                    params.R_cpe = R_cpe;
                    params.R_interue = R_interue(:,:,1:numCPE_tot,1:numCPE_tot);
                    params.R_ue = []; 
                    %stationary CPEs: no mobility layer; the fading draw
                    %count matches the cellular sample count (layers 2 x 3)
                    nbrFWADraws = nbrOfSnapshots*nbrOfRealizations;
                    rate_dl = zeros(numCPE_tot,nbrFWADraws);
                    for n = 1:nbrFWADraws
                        [~, ~, channel_dl_FWA, channel_est_dl_FWA] = computePhysicalChannels_sub6_MIMO(params);
                        %repeater-free MU-MIMO on the shared FWA band
                        rate_dl(:,n) = compute_link_rates_MIMO_mmse(params, channel_dl_FWA, channel_est_dl_FWA);
                    end
                    mean_rate_dl_FWA = mean(rate_dl,2);
                    if params.Band > 0 && num_rep_arr(idxnumrep) == 0
                        %raw FWA packing matrix (skipped if the cell phase
                        %consumed the whole band; first pool size only, with
                        %the activity operating point in the tag)
                        writematrix(rate_dl, strcat(packFolder,'/packing_FWA_',num2str(numCPE_all),'CPE_',num2str(activity_arr(idxActivity)),'act_',num2str(params.num_repeater_tot),'rep_',aID,'.csv'));
                    end
                    Band_FWA = params.Band;
                    %% FULL-BAND PER-TERMINAL SE COMPARISON (its own table):
                    %the gross physical-layer contrast - ALL cellular UEs
                    %thrown into the whole band vs ALL FWA CPEs thrown into
                    %the whole band. No orthogonal split, no demand or
                    %feasibility filtering, every terminal counted. The
                    %cellular per-UE SE is band-invariant under the
                    %constant-PSD convention, so the allocation-phase value
                    %is reused; the FWA side is recomputed explicitly on the
                    %full band (the num_rep = 0 fallow band can be zero,
                    %leaving nothing to rescale). Extra RNG draws here are
                    %safe: every num_rep iteration replays fadingState.
                    if num_rep_arr(idxnumrep) == 0
                        params_se = params; params_se.Band = Band;
                        rate_fb = zeros(numCPE_tot, nbrFWADraws);
                        for n = 1:nbrFWADraws
                            [~, ~, ch_fwa, ch_est_fwa] = computePhysicalChannels_sub6_MIMO(params_se);
                            rate_fb(:,n) = compute_link_rates_MIMO_mmse(params_se, ch_fwa, ch_est_fwa);
                        end
                        se_fwa_fb = sum(mean(rate_fb,2))/(Band*numCPE_tot);
                        se_cell_fb = sum(mean_rate_dl_cell)/(Band - Band_FWA)/(M_sectors*numUE);
                        seDir = fullfile('resultData','FWA_SE_comparison_logdet');
                        if not(isfolder(seDir)), mkdir(seDir), end
                        seFile = fullfile(seDir, strcat('se_comp_', ...
                            num2str(take_rate_arr(idxnumCPE)), '_', ...
                            num2str(activity_arr(idxActivity)), '_', aID, '.csv'));
                        fid_se = fopen(seFile,'w');
                        fprintf(fid_se, 'numCPE,numUE,se_cell_ue,se_fwa_cpe,multiple\n');
                        fprintf(fid_se, '%d,%d,%f,%f,%f\n', numCPE_tot, M_sectors*numUE, ...
                            se_cell_fb, se_fwa_fb, se_fwa_fb/se_cell_fb);
                        fclose(fid_se);
                    end
                    for idxrmin = 1:size(fwa_demand_configs,1)
                        %per-CPE absolute minimum rates for this config;
                        %the HOME rate labels the config in filenames/CSV
                        params.fwa_demand_rates = fwa_demand_configs(idxrmin, params.cpe_zone).';
                        params.fwa_demand_profile = fwa_demand_names{idxrmin}; %tags filenames/CSV
                        params.Band = Band_FWA;
                        %a CPE is SERVED iff its delivered rate meets its
                        %zone's advertised minimum for this demand config;
                        %the SE numerator sums the SERVED CPEs only, so
                        %numerator and denominator describe the same
                        %population (unserved residual rates would inflate
                        %the per-served-CPE SE, badly so when Band_FWA is
                        %small at num_rep = 0)
                        servedFWA = mean_rate_dl_FWA(:) >= params.fwa_demand_rates(:);
                        K_FWA_init = sum(servedFWA);
                        sum_FWA_rate = sum(mean_rate_dl_FWA(servedFWA));
                        params.Band = Band;
                        if (Band_FWA > 0) && (K_FWA_init > 0)
                            FWA_se_out = sum_FWA_rate/(Band_FWA*K_FWA_init);
                        else
                            FWA_se_out = 0; %no band or no served CPEs: report 0 instead of 0/0
                        end

                    %% Recording the Results (see the naming registry above):
                        %ONE file per seed, one row per sweep-grid point; the
                        %row prefix carries the CURRENT values of the swept
                        %parameters, in registry order
                        rateFolder = fullfile('resultData', naming.folder);
                        if not(isfolder(rateFolder))
                            mkdir(rateFolder)
                        end
                        resultFile = fullfile(rateFolder, naming.file);
                        writeHeader = ~isfile(resultFile);
                        currVals = {take_rate_arr(idxnumCPE), activity_arr(idxActivity), ...
                            params.num_repeater_tot, ...
                            params.fwa_demand_profile, params.repeat_gain, ...
                            params.SI_cancel_dB, params.r_min_cell/1e6, ...
                            params.loss_pc_cell, params.ISD, Band/1e6, params.HW_IMPAIRMENTS, ...
                            params.ncr_rank, params.num_repeater_per_cpe, ...
                            params.ncr_benefit_gate, params.fwa_pilot_contam, params.gI_site_aware, ...
                            logdet_generation};
                        fileID = fopen(resultFile,'a');
                        if writeHeader
                            fprintf(fileID, naming.header);
                        end
                        %cell_se_ue: cellular spectral efficiency PER UE (sum SE
                        %over the cellular band divided by the UE count) - the
                        %per-terminal counterpart of the per-CPE FWA_se column
                        fprintf(fileID, '%s,%d,%d,%d,%f,%f,%f\n', naming.rowfmt(currVals), ...
                            numCPE_all, M_sectors*numUE, K_FWA_init, Band_FWA, ...
                            sum(mean_rate_dl_cell)/(Band - Band_FWA)/(M_sectors*numUE), FWA_se_out);
                        fclose(fileID);
                    end
        end
    end
    end %idxnumCPE (FWA take-rate sweep)
tEnd = toc(tStart);
fprintf('Total runtime: %f seconds\n',tEnd)

