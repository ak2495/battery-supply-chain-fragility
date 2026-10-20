# Battery supply chain fragility
[![DOI](https://zenodo.org/badge/1427212137.svg)](https://doi.org/10.5281/zenodo.23503973)

This repository accompanies the paper:

<p align="center">
  <em>The Global Battery Supply Chain for the Green Transition Has Become More Fragile: Firm-Level Evidence from the Past Decade</em>
</p>

Amir H. Keshavarzzadeh, Qianyi Sun, Giacomo Vaccario, Tobias Hütner-Reisch, Kamiar Mohaddes, Dabo Guan and Jonathan Cullen.

The study examines changes in the structure and fragility of the global battery supply chain over 2015–2025 using firm-level supply relationships and simulated disruptions.

The code provides a simplified implementation of the NAICS-grouped Leontief propagation model used in the study. It demonstrates how distress spreads through supplier–buyer relationships and how firm-level susceptibility and transmission are calculated from individual firm failures. Supplier distress is averaged within each input category using transaction shares as weights. The most disrupted category determines the buyer’s input constraint.

The repository includes the propagation model, a script for calculating the impact matrix and firm-level measures, and a synthetic network with example results and numerical tests. It is provided to support understanding and reuse of the methodology. The empirical data and complete analysis underlying the paper’s figures and tables are not included.

## Files

| File | Purpose |
| --- | --- |
| `R/leontief.R` | Operator, propagation, impact matrix, susceptibility and transmission |
| `run_example.R` | Read a GraphML file and save the results |
| `data/` | Input requirements and a small synthetic graph |
| `tests/test_leontief.R` | Numerical and ordering tests, with analytic reference values |
| `output/` | Included synthetic example results and default destination for new results |
| `VALIDATION.md` | Validation scope, supplied run records and checks to run |
| `LICENSE` | MIT licence |

## Installation and use

Install R and the `igraph` and `Matrix` packages. No GPU or special hardware is needed for the synthetic example.

From R:

```r
install.packages(c("igraph", "Matrix"), repos = "https://cloud.r-project.org")
```

Open a terminal in the repository's main directory and run the tests and synthetic example:

```bash
Rscript tests/test_leontief.R
Rscript run_example.R
```

To analyse another network, pass the input file and output directory:

```bash
Rscript run_example.R /path/to/network.graphml output/sample
```

The default command reads `data/network.graphml` and writes to `output/`. Both scripts must be run from the repository root. Relative input and output paths are resolved from that working directory; absolute paths are also accepted.

The included graph is synthetic and does not represent the empirical battery supply chain. The proprietary data used for the study are not included.

The repository root contains `LICENSE`, `README.md`, `run_example.R` and `VALIDATION.md`, alongside the four directories `R/`, `data/`, `tests/` and `output/`. The input graph, its illustration and input documentation are in `data/`.

## Input graph

Edges point from suppliers to buyers. The graph must be directed, without parallel edges or self-links, and contain at least two firms. Each firm must have a numeric-digit NAICS code at the same digit level. The code checks this format, but does not check codes against an official classification. Each edge must have a finite, positive numeric transaction weight.

The default attributes are:

| Attribute | Meaning |
| --- | --- |
| Vertex `naics` | Industry classification, stored as a string |
| Edge `weight` | Transaction weight associated with the supplier–buyer relationship |
| Vertex `name` | Unique firm identifier; GraphML node `id` is used if `name` is absent |

Change the attribute settings at the top of `run_example.R` if the input file uses different names. The code checks the graph structure but does not clean or reclassify it. The NAICS edition, digit level, and the definition, units and source of the transaction weights should be documented with any supplied study sample.

## Method

### Input constraints

For each buyer and supplier NAICS category, transaction weights are divided by their total within that group. These shares are stored in a sparse matrix `P`, with one row per buyer–NAICS group and one column per supplier.

The function `leontief_operator()` calculates `P %*% h`, where `h` is the vector of firm distress, and then takes the maximum across each buyer's groups. A firm without suppliers receives zero input constraint.

### Distress propagation

The function `propagate_distress()` starts from the supplied shock vector and repeatedly applies:

```r
h_next <- pmax(h, leontief_operator(network, h))
```

Distress therefore cannot decrease during propagation. By default, iteration stops when `max(abs(h_next - h)) < 0.001`, or after 1,000 updates. The final iterate is returned in both cases, and a convergence flag records whether the change criterion was met. This criterion does not provide a bound on the distance to the exact limiting state.

### Impact matrix, susceptibility and transmission

A unit shock is applied to each firm in turn. The entry `M[i, j]` records the final distress at firm `i` following the failure of firm `j`. Rows and columns follow the same firm order, and the diagonal is one.

A firm's Leontief susceptibility is its row sum excluding the diagonal, divided by `n - 1`. It therefore measures the firm's average distress across the individual failures of all other firms. Every other firm is included with equal weight; no shock sources are sampled and no firm-size weights are applied.

A firm's distress transmission is its column sum excluding the diagonal, divided by `n - 1`. It measures the average distress that the firm's failure causes at all other firms. Susceptibility measures distress received; transmission measures distress caused. Both are dimensionless and lie between zero and one. These outputs are not monetary exposure estimates.

The impact matrix contains separate single-firm scenarios. Because the propagation operator is nonlinear, multiplying this matrix by a general shock vector does not, in general, reproduce the result of that shock. Use `propagate_distress()` directly for simultaneous or partial shocks.

For example, to apply a partial shock to firm A in the synthetic network:

```r
source("R/leontief.R")
g <- igraph::read_graph("data/network.graphml", format = "graphml")
network <- prepare_network(g)

shock <- numeric(network$n)
shock[match("A", network$ids)] <- 0.5

result <- propagate_distress(network, shock)
result$distress
```

Shock vectors must follow `network$ids` in order. Names on a supplied vector do not reorder its entries.

## Results

| Output | Contents |
| --- | --- |
| `impact_matrix.rds` | Sparse impact matrix with row and column firm identifiers, for R |
| `impact_matrix.mtx` | Sparse impact matrix in Matrix Market format |
| `matrix_firms.csv` | One-based matrix indices and firm IDs, shared by rows and columns |
| `susceptibility.csv` | Firm ID, NAICS code and susceptibility |
| `transmission.csv` | Firm ID and transmission |
| `convergence.csv` | Shocked firm, convergence flag, number of updates and final change |
| `run_settings.csv` | Input path and checksum, network size, parameters and elapsed calculation time |
| `session_info.txt` | R, operating system and package versions |

A warning is printed if any scenario reaches the iteration limit without satisfying the stopping criterion. Its final iterate is still saved. Repeating a run in the same output directory replaces the previous files.

To load the impact matrix in R:

```r
M <- readRDS("output/impact_matrix.rds")
```

For other languages, use `impact_matrix.mtx` together with `matrix_firms.csv`, because Matrix Market does not preserve firm identifiers. Rows identify receiving firms and columns identify shocked firms.

## Computational implementation

Both `P` and the impact matrix use compressed sparse column storage (`dgCMatrix`). The matrix `P` contains one stored value per edge, and no dense adjacency matrix is constructed. Impact matrix columns are assembled from the positive entries of each scenario without allocating a dense `n × n` matrix. No positive values are removed by a numerical cutoff. Distress and group constraints are stored as ordinary vectors.

A sparse supply network does not necessarily produce a sparse impact matrix. If most firms affect most other firms, storage still grows in proportion to `n²`, and sparse storage can require more memory than dense storage. Column assembly also requires temporary memory. This simplified implementation runs all shocks serially and stores the full result in memory.

## Synthetic example and validation

The synthetic graph contains eight firms, A–H. A and B share a NAICS category and supply D with weights 3 and 1. C supplies D from another category. D supplies E; E supplies F; F and G supply each other. E and G belong to different categories. H is isolated.

Failure of A gives D an input constraint of 0.75, which propagates to E, F and G. Failure of B gives 0.25 at the same firms, while failure of C gives 1. These values determine the analytic reference matrix embedded in `tests/test_leontief.R`; no separate `expected/` directory is needed.

Tests also cover a cycle with gradual convergence, the iteration limit, an edgeless graph, partial shocks, invalid weights, and reordering of edges and firms. The reference values were derived independently of the implementation. The supplied run environment and the scope of validation are documented in `VALIDATION.md`.

These tests assess the implementation on numerical examples. They do not establish agreement with empirical results generated by other scripts.

## Licence and citation

This simplified implementation is distributed under the MIT licence in `LICENSE`.

The accompanying manuscript title and author list are provided above. The code is available at https://github.com/ak2495/battery-supply-chain-fragility.

