# RQ1 Parameter Sweep

Exploratory LigandMPNN parameter tests carried out before settling on the
main reported pipeline (rq1/). Only one representative script per variant
type is kept here; near-duplicate values are noted below rather than
included individually, to avoid repetition.

- **Temperature sweep**: tested at 0.2, 0.3, 0.5, and 0.7 for both C2 and C3.
  Scripts for 0.2 and 0.7/0.5 (the tested extremes) are kept as
  representative examples; 0.3 (and 0.5 for C2) are omitted as they follow
  the identical script structure with only the --temperature value changed.
- **Seed sensitivity tests**: tested with seedtest_a through seedtest_f (C2)
  and seedtest_a through seedtest_c (C3), varying only the random seed.
  seedtest_a is kept as a representative example for each state.
- **Low-temperature runs**: tested with lowtemp_a and lowtemp_b (C2).
  lowtemp_a is kept as a representative example.
- **Mutation fix variants** (fixY75, fixY75L73): both kept, as these test
  genuinely different fixed-residue conditions rather than repeated
  parameter values (see Results 1.2, "Position 75" finding).
- **Docking-pose variants** (top1, top3, v2, base): all kept, as these
  represent distinct starting configurations rather than a parameter sweep.

None of the parameter sweep configurations outperformed the default
settings used in the main reported pipeline (rq1/); see dissertation
Methods for the final selected configuration.
