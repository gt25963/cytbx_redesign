# RQ1 Parameter Sweep

Exploratory LigandMPNN parameter tests carried out before settling on the main reported pipeline (rq1/). A representative script is kept per variant type; additional values/seeds mentioned below were explored on the HPC but not individually retained here. 

- **Temperature sweep**: `master_cytbx_C2_0.2.sh` / `master_cytbx_C3_0.2.sh` kept as representative examples. Other temperatures (0.1, 0.2, 0.3, 0.5, 0.7) were tested with the same script structure, varying only `--temperature`.
- **Seed sensitivity**: `master_cytbx_C2_0.5_seedtest_a.sh` / `master_cytbx_C3_0.3_seedtest_a.sh` kept as representative examples of the seed-sensitivity approach, varying only the random seed.
- **Mutation fix variant**: `master_cytbx_C2_0.5_fixY75.sh` kept as a representative example of fixing a specific residue during design.
- **Docking-pose variant**: `master_cytbx_C2_top1.sh` / `master_cytbx_C3_top1.sh` kept as representative examples; the main pipeline's default docking-pose selection is documented in rq1/.

Each `master_cytbx_*.sh` script has a matching `run_pipeline_*.sh` launcher. `run_4tool.sh` is a shared entry point used across these sweep variants.

>[!IMPORTANT]
>None of the parameter sweep configurations outperformed the default settings used in the main reported pipeline (rq1/); see report section *Methods* for the final selected configuration.
