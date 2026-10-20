# Synthetic example network

`network.graphml` is an artificial eight-firm supply chain network for demonstrating and testing the simplified code. `network_visualisation.png` illustrates the same graph. These files contain no empirical firm data.

A and B supply D from the same NAICS category with weights 3 and 1. C supplies D from another category with weight 2. D supplies E; E and G supply F from different categories; F supplies G. H is isolated. All remaining edge weights are 1. Weights are illustrative purchase values in arbitrary common units.

## Required structure

- Directed edges from suppliers to buyers, without parallel edges or self-links.
- At least two firms, with unique non-empty identifiers in vertex `name`, or GraphML node IDs if `name` is absent.
- Vertex `naics` containing digits only, with the same number of digits for all firms. Store codes as strings. Text labels and ranges are not accepted by the supplied validation.
- Edge `weight` containing a finite, positive numeric purchase value for each relationship. For A to B, this is the value buyer B purchases from supplier A.

## Configuration

Run `Rscript run_example.R` from the repository root to use this graph. To use another graph, run `Rscript run_example.R /path/to/network.graphml output/sample` from the same directory.

If attribute names differ, change `naics_attribute`, `weight_attribute` and, if needed, `id_attribute` in `run_example.R`. The identifier default uses `name` when present and otherwise `id`.

Changing the industry attribute selects a different field but does not translate classifications. A different coding scheme must satisfy the existing digits-only and equal-length checks, or those checks must be adapted explicitly. The grouping treats identical supplier codes as the same input category.

The code validates the input structure but does not clean data or verify classification codes against an official list. For an empirical application, document the classification edition and digit level and the definition, units and source of the weights.
