#!/usr/bin/env Rscript
# Run from the main project directory: Rscript tests/test_leontief.R
source("R/leontief.R")

check_equal <- function(x, y) {
  if (inherits(x, "Matrix")) x <- as.matrix(x)
  if (inherits(y, "Matrix")) y <- as.matrix(y)
  stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12)))
}
make_graph <- function(from, to, weight, codes, ids = names(codes)) {
  igraph::graph_from_data_frame(
    data.frame(from = from, to = to, weight = weight), directed = TRUE,
    vertices = data.frame(name = ids, naics = unname(codes))
  )
}

# Compare the full result with values worked out independently.
g <- igraph::read_graph("data/network.graphml", format = "graphml")
net <- prepare_network(g)
result <- compute_impact_matrix(net)
# Analytic reference: rows are receivers A-H, columns are failed firms A-H.
# A and B contribute 3/4 and 1/4 within their category at D; C is essential.
# These constraints pass through D -> E -> F <-> G. H is isolated.
expected <- rbind(
  c(1,    0,    0, 0, 0, 0, 0, 0),
  c(0,    1,    0, 0, 0, 0, 0, 0),
  c(0,    0,    1, 0, 0, 0, 0, 0),
  c(0.75, 0.25, 1, 1, 0, 0, 0, 0),
  c(0.75, 0.25, 1, 1, 1, 0, 0, 0),
  c(0.75, 0.25, 1, 1, 1, 1, 1, 0),
  c(0.75, 0.25, 1, 1, 1, 1, 1, 0),
  c(0,    0,    0, 0, 0, 0, 0, 1)
)
stopifnot(identical(net$ids, LETTERS[1:8]))
check_equal(result$M, expected)
check_equal(compute_susceptibility(result$M), c(0, 0, 0, 2, 3, 5, 5, 0) / 7)
check_equal(compute_transmission(result$M), c(3, 1, 4, 3, 2, 1, 1, 0) / 7)
stopifnot(all(result$diagnostics$converged), all(Matrix::diag(result$M) == 1))

# Within-category substitution and the bottleneck across categories.
h <- c(0.2, 0.6, 0.4, 0, 0, 0, 0, 0)
step <- leontief_operator(net, h)
check_equal(step[4L], 0.4)  # max(0.75 * 0.2 + 0.25 * 0.6, 0.4)
check_equal(leontief_operator(net, rep(0, 8)), rep(0, 8))
check_equal(leontief_operator(net, rep(1, 8))[c(1, 2, 3, 8)], rep(0, 4))

# This ordering exposed the numerator/denominator mismatch in the old code.
g_order <- make_graph(c("A", "B"), c("C", "C"), c(10, 20),
                      c(A = "331110", B = "325180", C = "335910"))
check_equal(leontief_operator(prepare_network(g_order), c(1, 0, 0)), c(0, 0, 1))

# Reordering edges and vertices must not change the result for any firm.
e <- igraph::as_data_frame(g, what = "edges")
v <- igraph::as_data_frame(g, what = "vertices")
g_reordered <- igraph::graph_from_data_frame(
  e[nrow(e):1L, ], directed = TRUE, vertices = v[nrow(v):1L, ]
)
r_reordered <- compute_impact_matrix(prepare_network(g_reordered))
check_equal(r_reordered$M[net$ids, net$ids], result$M)

# Distress stays in [0, 1] and cannot fall between successive updates.
h <- c(0.35, 0, 0.2, 0, 0, 0, 0, 0)
for (i in seq_len(20)) {
  next_h <- pmax(h, leontief_operator(net, h))
  stopifnot(all(next_h >= h), all(next_h >= 0 & next_h <= 1))
  h <- next_h
}

# A cycle can converge gradually when its suppliers share a category.
g_cycle <- make_graph(c("A", "C", "B"), c("B", "B", "C"), c(1, 1, 1),
                      c(A = "325180", B = "325180", C = "325180"))
cycle_net <- prepare_network(g_cycle)
cycle <- propagate_distress(cycle_net, c(1, 0, 0), tolerance = 1e-10)
stopifnot(cycle$converged, max(abs(cycle$distress - 1)) < 1e-8)
limited <- propagate_distress(cycle_net, c(1, 0, 0), max_iterations = 1L)
stopifnot(!limited$converged, limited$iterations == 1L)
check_equal(limited$distress, c(1, 0.5, 0))
check_equal(limited$final_change, 0.5)

# With no edges, the initial shock remains and susceptibility is zero.
g_empty <- igraph::make_empty_graph(2, directed = TRUE)
g_empty <- igraph::set_vertex_attr(g_empty, "name", value = c("A", "B"))
g_empty <- igraph::set_vertex_attr(g_empty, "naics", value = c("325180", "331110"))
empty <- compute_impact_matrix(prepare_network(g_empty))
check_equal(empty$M, diag(2))
check_equal(compute_susceptibility(empty$M), c(0, 0))

# The operator preserves order when the initial distress increases.
h_low <- seq(0, 0.7, length.out = 8)
h_high <- pmin(h_low + 0.2, 1)
stopifnot(all(leontief_operator(net, h_low) <= leontief_operator(net, h_high)))

# The input checks must reject invalid weights, without repairing them.
bad <- igraph::set_edge_attr(g, "weight", index = 1, value = -1)
stopifnot(inherits(try(prepare_network(bad), silent = TRUE), "try-error"))

# Verify sparse storage and group normalisation.
stopifnot(inherits(net$P, "dgCMatrix"), inherits(result$M, "dgCMatrix"))
check_equal(Matrix::rowSums(net$P), rep(1, nrow(net$P)))
stopifnot(length(net$P@x) == igraph::ecount(g))

# Compare sparse multiplication with an independent calculation over edges.
set.seed(123)
for (trial in seq_len(10)) {
  graph <- igraph::sample_gnp(15, 0.15, directed = TRUE, loops = FALSE)
  graph <- igraph::set_vertex_attr(graph, "name", value = as.character(seq_len(15)))
  codes <- sample(c("325180", "331110", "335910"), 15, replace = TRUE)
  graph <- igraph::set_vertex_attr(graph, "naics", value = codes)
  weights <- runif(igraph::ecount(graph), 1, 10)
  graph <- igraph::set_edge_attr(graph, "weight", value = weights)
  edge <- igraph::as_edgelist(graph, names = FALSE)
  h <- runif(15)
  direct <- numeric(15)
  for (buyer in seq_len(15)) {
    selected <- which(edge[, 2] == buyer)
    if (length(selected) > 0) {
      groups <- split(selected, codes[edge[selected, 1]])
      direct[buyer] <- max(vapply(groups, function(k) {
        sum(weights[k] * h[edge[k, 1]]) / sum(weights[k])
      }, numeric(1)))
    }
  }
  check_equal(leontief_operator(prepare_network(graph), h), direct)
}
# Check a partial shock through the full propagation routine.
partial <- propagate_distress(net, c(0.5, rep(0, 7)))
check_equal(partial$distress, c(0.5, 0, 0, rep(0.375, 4), 0))
stopifnot(partial$converged)

# Check exports without relying on or overwriting the bundled output files.
local({
  export_dir <- tempfile("leontief-export-")
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE))
  saveRDS(result$M, file.path(export_dir, "impact_matrix.rds"))
  invisible(Matrix::writeMM(result$M, file.path(export_dir, "impact_matrix.mtx")))
  restored <- readRDS(file.path(export_dir, "impact_matrix.rds"))
  check_equal(restored, expected)
  stopifnot(identical(rownames(restored), net$ids),
            identical(colnames(restored), net$ids))
  check_equal(Matrix::readMM(file.path(export_dir, "impact_matrix.mtx")), expected)
})
cat("All Leontief tests passed.\n")
