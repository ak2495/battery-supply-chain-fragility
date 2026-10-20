#!/usr/bin/env Rscript
# Run from the main project directory: Rscript run_example.R

if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("Install igraph first: install.packages('igraph')")
}
if (!requireNamespace("Matrix", quietly = TRUE)) {
  stop("Install Matrix first: install.packages('Matrix')")
}
source("R/leontief.R")

# The default input is the included synthetic example. Paths are relative to
# the working directory; run this script from the repository root.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) stop("Usage: Rscript run_example.R [graph_file] [output_dir]")
graph_file <- if (length(args) >= 1L) args[1L] else "data/network.graphml"
output_dir <- if (length(args) >= 2L) args[2L] else "output"

# These are the attribute names in the GraphML file.
naics_attribute <- "naics"
weight_attribute <- "weight"
id_attribute <- NULL  # Use 'name' when present, otherwise the GraphML 'id'.

# Numerical settings used in the methodology.
tolerance <- 0.001
max_iterations <- 1000L

if (!file.exists(graph_file)) {
  stop("Input file not found: ", graph_file,
       ". Run from the repository root or pass an existing GraphML path.")
}
g <- igraph::read_graph(graph_file, format = "graphml")
network <- prepare_network(g, naics_attribute, weight_attribute, id_attribute)

start <- proc.time()[["elapsed"]]
result <- compute_impact_matrix(network, tolerance, max_iterations)
elapsed <- proc.time()[["elapsed"]] - start
susceptibility <- compute_susceptibility(result$M)

# Rows are receiving firms. Columns are the firms whose failure is simulated.
# Write identifiers explicitly so the matrix order is preserved outside R.
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) stop("Could not create output directory: ", output_dir)
transmission <- compute_transmission(result$M)
write.csv(
  data.frame(
    firm_id = network$ids,
    transmission = as.numeric(transmission)
  ),
  file.path(output_dir, "transmission.csv"),
  row.names = FALSE
)
saveRDS(result$M, file.path(output_dir, "impact_matrix.rds"))
invisible(Matrix::writeMM(result$M, file.path(output_dir, "impact_matrix.mtx")))
write.csv(data.frame(matrix_index = seq_len(network$n), firm_id = network$ids),
          file.path(output_dir, "matrix_firms.csv"), row.names = FALSE)
write.csv(data.frame(firm_id = network$ids, naics = network$naics,
                     susceptibility = as.numeric(susceptibility)),
          file.path(output_dir, "susceptibility.csv"), row.names = FALSE)
write.csv(result$diagnostics, file.path(output_dir, "convergence.csv"), row.names = FALSE)
write.csv(data.frame(
  input_file = graph_file, input_md5 = unname(tools::md5sum(graph_file)),
  firms = network$n, edges = igraph::ecount(g),
  naics_attribute = naics_attribute, weight_attribute = weight_attribute,
  naics_digits = nchar(network$naics[1L]), tolerance = tolerance,
  max_iterations = max_iterations, converged_scenarios = sum(result$diagnostics$converged),
  elapsed_seconds = elapsed
), file.path(output_dir, "run_settings.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output_dir, "session_info.txt"))
cat(sprintf("Finished %d scenarios in %.3f seconds. Results: %s\n",
            network$n, elapsed, output_dir))
