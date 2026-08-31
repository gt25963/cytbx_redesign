# CytbX Redesign: Homo-Oligomeric Assembly and Cofactor Diversification

Computational pipeline for the iterative redesign of the de novo membrane cytochrome CytbX, developed as part of an MSc Bioinformatics dissertation (University of Bristol, 2026). Two parallel research objectives:

- **RQ1**: Controlled homo-oligomeric assembly (C2/C3) compatible with hexagonal 2D membrane lattices
- **RQ2**: Cofactor diversification via flavin mononucleotide (FMN) and ubiquinone-8 (Q8) binding

All computation was performed on the Isambard-AI HPC facility. Full methodological detail is in the accompanying dissertation; this README maps the codebase to that document.

## Repository Structure

- `rq1/` - Production pipeline scripts for RQ1 (oligomeric assembly)
  - `prescreening/` - Symmetry-constrained docking and combined-score prescreening (Table 1)
  - `compile/` - Score compilation across cycles
  - `master_cytbx_4tool*.sh` - Main iterative design drivers, one per lineage/cycle (Table 2)
- `rq2/` - Production pipeline scripts for RQ2 (cofactor diversification)
  - `fmn/` - FMN track drivers
  - `q8/` - Q8 track drivers
  - `u10/` - Discontinued ubiquinone-10 track (superseded by Q8; retained for transparency)
  - `params/` - Rosetta ligand parameter files (FMN.params, Q8.params, U10.params) and generation scripts
- `scripts/` - Core reusable utilities called throughout both pipelines
  - `scoring/` - Score extraction, compilation, and tool-input preparation (Chai-1, AF3, Boltz-2, ESM3)
  - `structure/` - Structural manipulation (chain splitting, oligomer building, distance measurement)
  - `verification/` - Direct coordinate/burial/clash verification (Discussion 4.3, 5.2, Table 6)
  - `prep/` - Ligand parameter generation
  - `figures/` - Plotting scripts for main-text figures. One representative script is kept per figure type rather than one script per individual figure number, since the underlying plotting logic is structurally identical across figures and differs only in the input data (see note below).
- `array_jobs/` - Individual SLURM submission scripts for each tool (RPXDock, LigandMPNN, Boltz-2, ESM3, Chai-1, AF3). Use these to resubmit a single pipeline step if a full run fails partway through, without rerunning the entire cycle.
- `exploratory/` - Trialled but not adopted in the reported pipeline
  - `rq1_parameter_sweep/` - LigandMPNN temperature sweep (0.2-0.7) and seed-sensitivity tests referenced in Methods but not individually detailed
  - `rosetta_alternative_approaches/` - FastRelax/FastDesign-based design routes (Methods: "Alternative approaches, control trials"), concluded as fair negative trials
  - `rq2_repose_stability/` - Early cofactor-placement and stability-scoring approaches superseded by the final LigandMPNN-based redesign pipeline
- `debug_oneoff/` - One-off diagnostic/recovery scripts used during development; not part of the core reported pipeline but retained for full transparency
- `example_outputs/` - Representative compiled score CSVs (one RQ1 C3 cycle, one RQ2 FMN cycle) showing the data format produced by `scripts/scoring/`. Full raw structure files (PDB/CIF, confidence JSONs) are not included due to size, but are available on request.

### Note on `scripts/figures/`

Where multiple figures share near-identical plotting logic (e.g. the RQ1 trajectory line plots, the RQ2 candidate-scatter comparisons), one representative script is kept per figure *type* rather than duplicated per figure number - the plotting code is the same, only the source CSVs differ. Current contents:

- `figure2B.py` / `figure2C.py` / `figure2D.py` - prescreening bar chart, per-cycle score distribution box-plots, iterative trajectory line plot (shown for the C3 track; C2 and the RQ2 tracks use the same plotting logic against different source CSVs)
- `figure4.py` - Chai-1 corrected vs aggregate score (covers both C2 and C3 panels in one script)
- `figure6.py` - Chai-1 vs AlphaFold 3 candidate scatter (covers all four tracks: C2, C3, FMN, Q8, in one script)
- `figure7.py` - FMN placement bias bar chart (design intent vs AF3/Boltz-predicted structure)
- `cofactor_panels.py` - FMN/Q8 chemical structure renders

Figures 1, 3, 9, 10, and the structural-panel components of Figures 5 and 8 were generated in PyMOL rather than matplotlib and are not scripted here. The data-loading logic for RQ2 figures (FMN/Q8) differs more substantially from the RQ1 examples than the plotting logic does - each cycle's source CSV has a different column-naming convention and, in some cases, requires lineage filtering - so the RQ1 scripts above are representative of the plotting approach but not a drop-in template for RQ2's data wiring.

## Key Methodological Notes

- **Chai-1 scores**: always extracted from `per_chain_pair_iptm[0][1]` in raw `.npz` output, not the summary `iptm` field (which is inflated by cofactor prediction confidence)
- **Cofactor placement bias**: AF3 and Boltz-2 systematically relocate FMN/Q8 to the incorrect haem site relative to design intent (Discussion 5.1). All coordination distances reported in the dissertation were measured directly on packed structures, not predicted structures - see `scripts/verification/`
- **Burial ratio correction**: FMN candidate selection metrics were found to reflect AF3's mispredicted structure in earlier cycles; corrected packed-structure verification is in `scripts/verification/burial_and_coordination_check.py` (Discussion 5.2, Table 6)

## Software Dependencies

LigandMPNN, PyRosetta, RPXDock, Boltz-2, ESM3, Chai-1, AlphaFold 3, RDKit, Biopython. See dissertation Methods for full environment/version details.

## Data Availability

Raw structure files (PDB/CIF outputs from each design cycle) are not included in this repository due to size, but are available on request.
