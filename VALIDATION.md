# Validation

This repository provides a simplified implementation and a synthetic numerical example. The checks below assess the implementation; they do not establish equivalence to the full empirical analysis or reproduce the manuscript's figures and tables.

## Supplied example results

The files in `output/` were supplied with this repository. They correspond to the eight-firm, seven-edge graph in `data/network.graphml`, using tolerance `0.001` and at most 1,000 updates per scenario.

The supplied `output/session_info.txt` records:

- R 4.4.1
- igraph 2.2.2
- Matrix 1.7-0
- macOS 26.5.2, aarch64-apple-darwin20

The supplied `output/run_settings.csv` records an impact-matrix calculation time of approximately 0.033 seconds. This excludes graph loading and file export and is not a performance benchmark. All eight scenarios are marked as converged in `output/convergence.csv`, with final change zero.

The original run records are retained unchanged. They describe the supplied run, not a new execution of the revised scripts.

## Checks performed during this review

The GraphML input checksum, firm ordering, Matrix Market impact matrix, susceptibility, transmission and convergence records were independently checked against an edge-based calculation and the analytic eight-firm example. The supplied RDS matrix was also checked against the Matrix Market matrix using a Python RDS reader.

R was unavailable in the review environment, so the revised R scripts and test suite have not been executed here. Static review and independent checks of the supplied output do not replace running the commands below before release.

## Running the tests

From the repository root, with `igraph` and `Matrix` installed:

```bash
Rscript tests/test_leontief.R
Rscript run_example.R
```

A successful test run prints `All Leontief tests passed.` The example command reads `data/network.graphml` and overwrites the result files in `output/`, including fresh run settings and session information. To retain the supplied outputs, use a separate destination:

```bash
Rscript run_example.R data/network.graphml output/check
```

Analytic reference values are embedded in `tests/test_leontief.R`, so tests do not depend on a missing reference directory or on previously generated outputs. The suite checks the impact matrix, susceptibility, transmission, within-category weighted constraints, bottlenecks across categories, partial shocks, isolation, an edgeless graph, gradual convergence through a cycle, the iteration limit, persistence, bounds, order preservation and invariance to edge and vertex ordering.

It also checks sparse storage, group normalisation, invalid weights, comparison with an independent edge-based operator on ten random graphs, and RDS and Matrix Market round trips using temporary files.

The stopping criterion limits the change between successive updates; it does not bound the error relative to the exact limiting state. Non-converged scenarios retain their final iterate and produce a warning.
