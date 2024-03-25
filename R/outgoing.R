#' Creates an outgoing tree from a given undirected treee.
#' @param tree Input tree, undirected.
#' @param dir_tree Directionalities that must be present.
#' @param a The node being inspected.
#' @param b The neighbor being inspected.

outgoing = function(tree, dir_tree = NULL, a = NULL, b = 1) {
	p = vcount(tree)
	if (is.null(dir_tree)) {
		dir_tree = make_empty_graph(p)
		nn = neighbors(tree, 1) 
		if (length(nn) == 0) return(dir_tree)
		for (j in nn) {
			dir_tree = add_edges(dir_tree, c(1, j))
		}
		for (j in nn) {
			dir_tree = outgoing(tree, dir_tree, a = 1, b = j)
		}
	}
	else {
		nn = setdiff(neighbors(tree, b), a)
		if (length(nn) == 0) return(dir_tree)
		for (j in nn) {
			dir_tree = add_edges(dir_tree, c(b, j))
		}
		for (j in nn) {
			dir_tree = outgoing(tree, dir_tree, a = b, b = j)
		}		
	}
	return(dir_tree)
}

