# Simplified NAICS-grouped Leontief propagation accompanying the manuscript.
#
# Each node represents a firm, and each edge points from a supplier to a buyer.
# Distress ranges from 0 (unaffected) to 1 (complete failure).
#
# For each buyer, supplier distress is averaged within each industry,
# weighted by purchase value. The largest industry-level constraint
# determines the buyer's overall input constraint.
#
# All input vectors must follow the graph's vertex order.


# Check the graph and prepare the supplier groups and purchase shares.
prepare_network <- function(g, naics_attribute = "naics",
                            weight_attribute = "weight", id_attribute = NULL) {

  # Reject undirected graphs, duplicate edges and links from a firm to itself (loop).
  if (!igraph::is_directed(g) || !igraph::is_simple(g)) {
    stop("Use a directed graph without parallel edges or self-links.")
  }
  n <- igraph::vcount(g)
  if (n < 2L) stop("Use at least two firms to calculate susceptibility.")

  # Unless an identifier attribute is specified, use "name" if available;
  # otherwise, try "id".
  if (is.null(id_attribute)) {
    id_attribute <- if ("name" %in% igraph::vertex_attr_names(g)) "name" else "id"
  }

  # Read firm identifiers, industry codes and purchase values.
  ids <- as.character(igraph::vertex_attr(g, id_attribute))
  naics <- as.character(igraph::vertex_attr(g, naics_attribute))
  w <- igraph::edge_attr(g, weight_attribute)

  # Every firm must have a non-empty identifier that is used only once.
  if (length(ids) != n || anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Firm identifiers must be present, non-empty and unique.")
  }

  # Require industry codes containing digits only, all at the same level
  # of detail, for example six digits for every firm.
  if (length(naics) != n || anyNA(naics) || any(!grepl("^[0-9]+$", naics))) {
    stop("Each firm must have a NAICS code containing digits only.")
  }
  if (length(unique(nchar(naics))) != 1L) {
    stop("Use the same NAICS digit level for all firms.")
  }

  # Every existing edge must have a finite purchase value greater than zero.
  if (igraph::ecount(g) > 0L &&
      (!is.numeric(w) || length(w) != igraph::ecount(g) ||
       any(!is.finite(w)) || any(w <= 0))) {
    stop("Each edge must have a finite, positive numeric weight.")
  }

  # Extract the supplier and buyer vertex positions for each edge.
  edges <- igraph::as_edgelist(g, names = FALSE)
  supplier <- edges[, 1L]
  buyer <- edges[, 2L]

  # Group each buyer's suppliers by their industry code.
  # Suppliers belong to the same group only if they supply the same buyer
  # and have the same NAICS code.
  key <- paste(buyer, naics[supplier], sep = ":")
  group <- match(key, unique(key))

  # Record which buyer each group belongs to, in group-number order.
  group_buyer <- buyer[!duplicated(group)]

  # Divide each purchase value by the total purchases in its group.
  # For example, purchases of 30 and 70 from two suppliers in the same
  # group give shares of 0.3 and 0.7.
  share <- numeric(length(w))
  if (length(w) > 0L) {
    total <- as.numeric(rowsum(w, group, reorder = TRUE))
    if (any(!is.finite(total))) stop("Group weights are too large to sum.")
    share <- w / total[group]
  }

  # Store the shares in matrix P: rows are buyer-industry groups,
  # columns are suppliers, and entries are purchase shares.
  P <- Matrix::sparseMatrix(
    i = group, j = supplier, x = share,
    dims = c(length(group_buyer), n), repr = "C"
  )

  list(n = n, ids = ids, naics = naics, P = P, group_buyer = group_buyer)
}


# Calculate each buyer's input constraint from current supplier distress.
leontief_operator <- function(network, h) {

  # Require one distress value between 0 and 1 for every firm.
  stopifnot(is.numeric(h), length(h) == network$n,
            all(is.finite(h)), all(h >= 0 & h <= 1))

  # Start with zero constraints, including for firms without suppliers.
  constraint <- numeric(network$n)

  if (nrow(network$P) > 0L) {

    # Calculate average supplier distress within each buyer-industry group,
    # using purchase shares as weights.
    group_constraint <- as.numeric(network$P %*% h)

    # A buyer is constrained by its most disrupted industry group.
    # Inputs from other industries cannot compensate for that constraint.
    buyer_constraint <- tapply(group_constraint, network$group_buyer, max)

    # Place each result at the corresponding buyer's vertex position.
    constraint[as.integer(names(buyer_constraint))] <- as.numeric(buyer_constraint)
  }

  # Prevent rounding errors from producing values above 1,
  # and label the results with firm identifiers.
  setNames(pmin(constraint, 1), network$ids)
}


# Repeatedly pass distress from suppliers to buyers until changes are small
# enough, or the maximum number of iterations is reached.
propagate_distress <- function(network, shock, tolerance = 0.001,
                              max_iterations = 1000L) {

  # shock gives the initial distress of every firm, between 0 and 1.
  # It may affect one firm or several firms.
  # tolerance is the stopping threshold for the largest change in distress.
  stopifnot(is.numeric(shock), length(shock) == network$n,
            all(is.finite(shock)), all(shock >= 0 & shock <= 1),
            length(tolerance) == 1L, is.finite(tolerance), tolerance > 0,
            length(max_iterations) == 1L, is.finite(max_iterations),
            max_iterations >= 1, max_iterations == floor(max_iterations))

  h <- setNames(as.numeric(shock), network$ids)
  converged <- FALSE

  for (iteration in seq_len(max_iterations)) {

    # Keep whichever is larger: existing distress or the new input constraint.
    h_next <- pmax(h, leontief_operator(network, h))

    # Measure the largest change across all firms in this iteration.
    change <- max(abs(h_next - h))
    h <- h_next

    # Stop when every firm's change is below the chosen tolerance.
    if (change < tolerance) {
      converged <- TRUE
      break
    }
  }

  # Return the final distress values and information about the calculation.
  # If the iteration limit was reached, converged remains FALSE.
  list(distress = h, converged = converged, iterations = iteration,
       final_change = change)
}


# Run a separate complete-failure scenario for every firm.
# In the resulting matrix M, entry [i, j] is the final distress of firm i
# when firm j is the initial source of failure.
compute_impact_matrix <- function(network, tolerance = 0.001,
                                  max_iterations = 1000L) {

  n <- network$n

  # Store affected firm positions and their distress values for each scenario.
  # Zero values are omitted to save memory; all positive values are retained.
  rows <- vector("list", n)
  values <- vector("list", n)

  # Record whether each scenario converged, how many iterations it took,
  # and the largest distress change in its final iteration.
  diagnostics <- data.frame(
    shock_firm = network$ids, converged = logical(n),
    iterations = integer(n), final_change = numeric(n)
  )

  for (j in seq_len(n)) {

    # Start a new scenario with firm j fully failed and all others unaffected.
    shock <- numeric(n)
    shock[j] <- 1

    result <- propagate_distress(network, shock, tolerance, max_iterations)

    # Save the positive distress values and convergence information.
    rows[[j]] <- which(result$distress > 0)
    values[[j]] <- unname(result$distress[rows[[j]]])
    diagnostics$converged[j] <- result$converged
    diagnostics$iterations[j] <- result$iterations
    diagnostics$final_change[j] <- result$final_change
  }

  # Assemble the scenarios as columns of a sparse matrix.
  # p marks where each column's entries begin in the stored values.
  M <- Matrix::sparseMatrix(
    i = unlist(rows, use.names = FALSE), p = c(0, cumsum(lengths(rows))),
    x = unlist(values, use.names = FALSE), dims = c(n, n),
    dimnames = list(network$ids, network$ids), repr = "C"
  )

  # Warn if any scenario stopped before meeting the convergence criterion.
  # Its final approximation is still included in M.
  if (any(!diagnostics$converged)) {
    warning(sum(!diagnostics$converged),
            " scenarios reached the iteration limit. Inspect the returned diagnostics.",
            call. = FALSE)
  }

  list(M = M, diagnostics = diagnostics)
}


# Calculate each firm's average distress across failures of all other firms.
compute_susceptibility <- function(M) {

  stopifnot(inherits(M, "sparseMatrix"), nrow(M) == ncol(M), nrow(M) > 1L)

  # Each row contains one firm's distress across all failure scenarios.
  # Remove its own-failure scenario, which lies on the diagonal,
  # then divide by the number of other firms.
  (Matrix::rowSums(M) - Matrix::diag(M)) / (nrow(M) - 1L)
}


# Calculate the average distress a firm's failure causes at all other firms.
compute_transmission <- function(M) {
  stopifnot(inherits(M, "sparseMatrix"), nrow(M) == ncol(M), nrow(M) > 1L)

  # Each column is one failure scenario; exclude the initially failed firm.
  (Matrix::colSums(M) - Matrix::diag(M)) / (ncol(M) - 1L)
}
