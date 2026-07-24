# Can the Revenue Potential for 6G Fixed Wireless Access Match Mobile Cellular?

Simulation and figure-generation code for the paper:

> **Can the Revenue Potential for 6G Fixed Wireless Access Match Mobile Cellular?**
> Soumyadeep Datta, Ish Kumar Jain, Rohit Budhiraja, Pei Liu, and Shivendra S. Panwar.
> *Submitted to IEEE Communications Magazine, 2026.*

**Authors and affiliations**
- Soumyadeep Datta — NYU Tandon School of Engineering, USA, and IIT Kanpur, India
- Ish Kumar Jain — Rensselaer Polytechnic Institute (RPI), USA
- Rohit Budhiraja — IIT Kanpur, India
- Pei Liu — NYU Tandon School of Engineering, USA
- Shivendra S. Panwar — NYU Tandon School of Engineering, USA

## Overview

Fixed Wireless Access (FWA) has become a major revenue source for 5G carriers, but its
revenue per gigabyte is roughly an order of magnitude below cellular. This work studies
whether FWA can be made competitive on the more decision-relevant metric of **revenue per
unit spectrum ($/Hz)** — a scaled product of revenue per GB and spectral efficiency (SE) —
by exploiting the larger form factor and stationary location of FWA customer premises
equipment (CPEs). The simulator models a combined downlink FWA–cellular deployment and
quantifies three levers: (i) higher-order MU-MIMO with inter-cell interference cancellation
for FWA CPEs; (ii) network-controlled repeaters (NCRs) collocated with FWA CPEs that raise
cellular SE and thereby free spectrum for FWA; and (iii) joint FWA–cellular spectrum
allocation.

The system model is a 3GPP TR 38.901 suburban-macro (SMa) deployment: two tri-sectored
macro sites over a zoned suburban layout (apartment complexes, single-family homes, a strip
mall, office parks, and a vegetation buffer), with a building-inventory CPE supply, a
busy-hour activity model, mobility-aware imperfect CSI (per-class sounding periodicity with
Jakes channel aging), and 3GPP-compliant hardware impairments.

## Repository layout

**Core simulation (MATLAB, tested on R2025b)**
- `SimulationMain.m` — entry point. Builds the drop, runs the cellular and FWA phases over
  the parameter sweep, and writes per-seed result tables and raw packing matrices.
- `generateSetup.m` — per-drop setup: SMa path loss / LOS / O2I, sectored antenna gains,
  spatial-correlation Cholesky factors, and the mobility state.
- `computePhysicalChannels_sub6_MIMO.m` — channel realizations with imperfect CSI (CSI-RS /
  SRS pilot overhead and per-class Jakes aging).
- `compute_link_rates_OFDM.m` / `compute_link_rates_OFDM_wi_repeater.m` — cellular (OFDMA)
  rates, without and with NCR assistance.
- `compute_link_rates_MIMO_mmse.m` / `compute_link_rates_MIMO_mmse_wi_repeater.m` — FWA
  MU-MIMO MMSE rates, without and with NCR-aided hybrid scheduling.
- `computeUtility.m` — served-CPE and utility computation against per-zone demand.
- `sweepNaming.m` — automated result-naming registry: swept parameters become table
  columns; constants name the output folder (`FWA_const_<names>`) and file
  (`results_<values>_<seed>.csv`).

**Post-processing (run offline, after all HPC array tasks finish)**
- `dataProcess/combineData.m` — aggregates per-seed CSVs into per-configuration summaries
  and reduces the packing matrices to safe-load-fraction curves.
- `dataProcess/plotData.m` — renders the paper figures in IEEE format.

**HPC**
- `submit.sbatch` — SLURM array entry point (one array task per random seed).

**Figure generators**
- `diagrams/make_deployment_3d.m` — deployment scene (2.5D icon map); also exports scene geometry.
- `diagrams/make_deployment_combined.m` — the system-model figure: deployment scene with the NCR and interference-nulling mechanisms drawn in.
- `diagrams/render_scene.py` — headless Blender render of the exported scene.
- `diagrams/make_mechanism_figs.m` — standalone NCR-support and interference-nulling vignettes.

## Running

Single local run:

```matlab
SimulationMain   % writes into resultData/FWA_const_*/
```

Multi-seed campaign on a SLURM cluster:

```bash
sbatch --array=1-N submit.sbatch                     # one seed per array task
matlab -batch "run('dataProcess/combineData.m')"     # after all tasks finish
matlab -batch "run('dataProcess/plotData.m')"        # render figures
```

Deployment and mechanism figures:

```bash
matlab -batch "run('diagrams/make_deployment_3d.m')"
/Applications/Blender.app/Contents/MacOS/Blender -b -P diagrams/render_scene.py   # optional
matlab -batch "run('diagrams/make_mechanism_figs.m')"
```

## Requirements

- MATLAB R2025b (Communications and Statistics toolboxes)
- SLURM for the multi-seed campaign
- Blender (optional) for the photorealistic deployment render

## Citation

This work is under review at IEEE Communications Magazine. Please check back for the final
citation and DOI. In the interim:

```bibtex
@article{datta2026fwa,
  author  = {Datta, Soumyadeep and Jain, Ish Kumar and Budhiraja, Rohit and Liu, Pei and Panwar, Shivendra S.},
  title   = {Can the Revenue Potential for 6G Fixed Wireless Access Match Mobile Cellular?},
  journal = {IEEE Communications Magazine (submitted)},
  year    = {2026}
}
```

## Acknowledgment

This work was supported in part by the "Next Generation Wireless Research and
Standardization on 5G and Beyond" project of the Ministry of Electronics and Information
Technology (MeitY), Govt. of India; the NY State Center for Advanced Technology in
Telecommunications (CATT); NYU Wireless; and the U.S. National Science Foundation (NSF)
through RINGS under Grant CNS-2148309.
