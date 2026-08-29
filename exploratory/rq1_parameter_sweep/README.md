# RQ1 Parameter Sweep

Exploratory LigandMPNN parameter tests carried out before settling on the
main reported pipeline (rq1/). One representative script per variant type
is kept; near-duplicate values and repeat configurations are noted below
rather than included individually, to avoid repetition.

- **Temperature sweep**: tested at 0.2, 0.3, 0.5, and 0.7 for C2, and 0.2,
  0.3, and 0.5 for C3. The 0.2 script is kept as a representative example
  for each state; the other values follow the identical script structure
  with only the --temperature value changed.
- **Seed sensitivity tests**: tested with seedtest_a through seedtest_f (C2)
  and seedtest_a through seedtest_c (C3), varying only the random seed.
  seedtest_a is kept as a representative example for each state.
- **Low-temperature runs**: tested with lowtemp_a and lowtemp_b (C2), both
  omitted here as they follow the same structure as the kept temperature
  sweep example.
- **Mutation fix variants**: tested fixY75 and fixY75L73 (C2). fixY75 is
  kept as a representative example; fixY75L73 follows the same approach
  with an additional fixed residue.
- **Docking-pose variants**: tested top1, top2 (main pipeline default), top3,
  v2, and the unlettered base configuration for both C2 and C3. top1 is
  kept as a representative example for each state.

None of the parameter sweep configurations outperformed the default
settings used in the main reported pipeline (rq1/); see dissertation
Methods for the final selected configuration.
