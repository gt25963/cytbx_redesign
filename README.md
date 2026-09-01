# CytbX Redesign: Homo-Oligomeric Assembly and Cofactor Diversification

Computational pipeline for the iterative redesign of the *de novo* membrane cytochrome CytbX, developed as part of an MSc Bioinformatics dissertation (University of Bristol, 2026). 
Two parallel research objectives:

- **RQ1**: Controlled homo-oligomeric assembly (C2/C3) compatible with hexagonal 2D membrane lattices
- **RQ2**: Cofactor diversification via flavin mononucleotide (FMN) and ubiquinone-8 (Q8) binding

All computation was performed on the Isambard-AI HPC facility. Full methodological detail is in the Methods section of the dissertation; this README maps the codebase to that document and examples of scripts used during the project.

## Repository Structure

- `rq1/` - Production pipeline RQ1 (oligomeric assembly)
  - `prescreening/` - RPXDock prescreening + rescoring (Table 1, Figure 2B)
  - `compile/` - Score compilation across cycles (C3)
  - `master_cytbx_4tool*.sh` - Main iterative design pipelines, one script per lineage/cycle (C2, C3, + several seed variants from later cycles) (Table 2)
- `rq2/` - Production pipeline RQ2 (cofactor diversification)
  - `fmn/` - FMN track 
  - `q8/` - Q8 track 
  - `u10/` - Discontinued ubiquinone-10 track (superseded by Q8; retained for transparency)
  - `params/` - Rosetta ligand parameter files (FMN.params, Q8.params, U10.params) and their generation/relaxation scripts
- `scripts/` - Reusable utilities called throughout both pipelines
  - `scoring/` - Score extraction, compilation, and tool-input preparation (Chai-1, AF3, Boltz-2, ESM3)
  - `structure/` - Structural manipulation (chain splitting, oligomer building, distance measurement)
  - `verification/` - Direct coordinate/burial verification (Discussion 5.2, Table 6)
  - `prep/` - Ligand parameter generation
  - `figures/` - Plotting scripts for dissertation figures (see note below)
  - `superseded/` - Earlier failed or not optimal versions (kept for transparency)
- `array_jobs/` - Individual SLURM submission scripts for each tool (RPXDock, LigandMPNN, Boltz-2, ESM3, Chai-1, AF3). Used these to resubmit a single pipeline step if a full run failed partway through, without rerunning the entire cycle.
- `exploratory/` - Trialled but not adopted in the reported pipeline
  - `rosetta_alternative_approaches/` - Concluded as fair trials but did not improve what was settled on in the end
    - `approach2/` - FastRelax + LigandMPNN
    - `approach3/` - RosettaMP FastDesign 
  - `rq1_parameter_sweep/` - LigandMPNN temperature sweep (0.2-0.7) and seed-sensitivity tests referenced in Methods (representative scripts - not detailed individually) 
  - `rq2_repose_stability/` - Early cofactor-placement + stability-scoring approaches superseded by the final LigandMPNN-based redesign pipeline
- `debug_oneoff/` - One-off diagnostic/recovery scripts used during development; not part of the core reported pipeline but retained for full transparency
- `example_outputs/` - Representative compiled score CSVs (one RQ1 C3 cycle, one RQ2 FMN cycle) showing the data format produced by `scripts/scoring/`.

>[!IMPORTANT]
>Full raw structure files (PDB/CIF, confidence JSONs) are not included due to size, but are available on request.

> [!NOTE]
> Ubiquinone-10 (U10/Q10) superseded the Q8 track. All scripts that contain U10 have been kept in for transparency, but know that for the final dissertation and project, Q8 is the correct cofactor.
 
### Note on `scripts/figures/`

Where multiple figures share near-identical plotting logic (e.g. the RQ1 trajectory line plots, the RQ2 candidate-scatter comparisons), one representative script is kept per figure type - the plotting code is the same, only the source CSVs differ:

- `figure2B.py` / `figure2C.py` / `figure2D.py` - prescreening bar chart, per-cycle score distribution box-plots, iterative trajectory line plot (shown for the C3 track; C2 and the RQ2 tracks use the same plotting logic against different source CSVs)
- `figure4.py` - Chai-1 corrected vs aggregate score (covers both C2 and C3 panels)
- `figure6.py` - Chai-1 vs AlphaFold 3 candidate scatter (covers all four tracks: C2, C3, FMN, Q8)
- `figure7.py` - FMN placement bias bar chart (design intent vs AF3/Boltz-predicted structure)
- `cofactor_panels.py` - FMN/Q8 chemical structure renders

Figures 1, 3, 9, 10, and the structural-panel components of Figures 5 and 8 were generated in PyMOL rather than matplotlib and are not scripted here. For RQ2 each cycle's source CSV has a different column-naming convention and, in some cases, requires lineage filtering - therefore the RQ1 scripts above are representative of the plotting approach but not a drop-in template.

## Key Methodological Notes

- **Chai-1 scores**: always extracted from *per_chain_pair_iptm[0][1]* in raw *.npz* output, not the summary *iptm* field (which is inflated by cofactor prediction confidence)
- **Cofactor placement bias**: AF3 and Boltz-2 systematically relocate FMN/Q8 to the incorrect haem site relative to design intent (Discussion 5.1). All coordination distances reported in my dissertation were measured directly on packed structures, not predicted structures - see *scripts/verification/*
- **Burial ratio correction**: FMN candidate selection metrics were found to reflect AF3's structure in earlier cycles; corrected packed-structure verification can be found in *scripts/verification/fmn_step7_compile_v2.py* (Discussion 5.2, Table 6)
- **Score normalisation**: *rq1/prescreening/rescore_simple.py* provides the raw un-normalised averaging used for all cross-state comparisons (Table 1)

## Software Dependencies

LigandMPNN, PyRosetta, RPXDock, Boltz-2, ESM3, Chai-1, AlphaFold 3, RDKit, Biopython. See dissertation Methods for full environment/version details.

## Data Availability
>[!IMPORTANT]
>Raw structure files (PDB/CIF/confidence JSONs/outputs from each design cycle) are not included in this repository due to size, but are available on request.
