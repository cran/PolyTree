#' Causal Polytree Estimation
#'
#' Estimates directed causal polytree from data, using algorithm developed in Chatterjee and Vidyasagar (2022).
#' @param x Data matrix, whose rows are i.i.d. data vectors generated from the model.
#' @return A directed polytree estimated from the input data, as an igraph object. 
#' @references Sourav Chatterjee and Mathukumalli Vidyasagar (2022). Estimating large causal polytrees from small samples. Available at <https://arxiv.org/abs/2209.07028>
#' @examples
#' p <- 10
#' n <- 200
#' x <- matrix(nrow = n, ncol = p)
#' for (i in 1:n) {
#'    x[i,1] = rnorm(1)
#'    for (j in 2:p) {
#'       x[i,j] = (x[i,j-1] + rnorm(1))/sqrt(2)
#'    }
#' }
#' p <- polytree(x)
#' @export
#' @import igraph
#' @import FOCI


polytree = function(x) {
	p = ncol(x)
	g = condeptree(x)
	dir_g = make_empty_graph(p, directed = TRUE) # Empty graph on p vertices
	# We will now implement steps 1 and 2 of the algorithm.
	determined = 1 # Number of edges whose directionalities 
	               # are determined in a step of the iteration
	iter = 0
	while (determined > 0) {
		determined = 0
		for (i in 1:p) {
			# i = the vertex being examined
			neighborhood = as.vector(neighbors(g, i)) # The set of neighbors of i in g
			dir_nbhd = as.vector(neighbors(dir_g, i, mode = "all")) # Set of all neighbors of i in the directed graph
			dir_nbhd_in = as.vector(neighbors(dir_g, i, mode = "in")) # Set of incoming neighbors of i
			if (length(dir_nbhd_in) == 0) {
				if (iter == 0) {
					if (length(neighborhood) > 0) {
						s = 0
						for (j in sort(neighborhood)) {
							for (k in sort(neighborhood)) {
								if (k > j) {
									xi1 = xicorln(x[,j], x[,k])   #codec(x[,k], x[,j])
									xi2 = codec(x[,k], x[,j], x[,i])
									if (xi2 >= xi1) {
										s = 1
										nn = as.vector(neighbors(dir_g, i, mode = "all"))
										if (is.element(j, nn) == FALSE) {
											dir_g = add_edges(dir_g, c(j,i))  # Incoming edge from j to i
											determined = determined + 1
										}
										if (is.element(k, nn) == FALSE) {
											dir_g = add_edges(dir_g, c(k,i))  # Incoming edge from k to i
											determined = determined + 1
										}
									}
								}
								if (s == 1) break
							}
							if (s == 1) break
						}
					}
				}
			}
			else {
				j = dir_nbhd_in[1] # This chooses a neighbor with an incoming edge to i
				unassigned = setdiff(neighborhood, dir_nbhd) # These are neighbors with 
															# unassigned directionalities of edges to i
				if (length(unassigned) > 0) {
					for (k in unassigned) {
						xi1 = xicorln(x[,j], x[,k])  #codec(x[,k], x[,j])
						xi2 = codec(x[,k], x[,j], x[,i])
						if (xi2 >= xi1) {
							dir_g = add_edges(dir_g, c(k,i))  # Incoming edge from k to i
							determined = determined + 1
						}
						else {
							dir_g = add_edges(dir_g, c(i,k))  # Outgoing edge from i to k
							determined = determined + 1
						}
					}
				}
			}
		}
		iter = 1 # Marks the end of the first iteration
	}
	# Next, we are going to implement step 3 of the algorithm.
	determined = 1 # Number of edges whose directionalities 
	               # are determined in a step of the iteration
	while (determined > 0) {
		determined = 0
		for (i in 1:p) {
			# i = the vertex being examined
			neighborhood = as.vector(neighbors(g, i)) # The set of neighbors of i in g
			dir_nbhd = as.vector(neighbors(dir_g, i, mode = "all")) # Set of all neighbors of i in the directed graph
			dir_nbhd_in = as.vector(neighbors(dir_g, i, mode = "in")) # Set of incoming neighbors of i			
			unassigned = setdiff(neighborhood, dir_nbhd) # These are neighbors with 
															# unassigned directionalities of edges to i
			if (length(dir_nbhd_in) > 0) {
				if (length(unassigned) > 0) {
					for (k in unassigned) {
						dir_g = add_edges(dir_g, c(i,k))  # Outgoing edge from i to k
						determined = determined + 1
					}
				}				
			}
		}
	}
	# Create "outgoing tree", i.e., tree having no vertex with more than one incoming edge
	outgoing_tree = outgoing(g)
	for (i in 1:p) {
		nn_in = setdiff(neighbors(outgoing_tree, i, mode = "in"), neighbors(dir_g, i, mode = "all")) # Incoming neighbors not already in dir_g
		nn_out = setdiff(neighbors(outgoing_tree, i, mode = "out"), neighbors(dir_g, i, mode = "all")) # Outgoing neighbors not already in dir_g
		if (length(nn_in) > 0) {
			for (j in nn_in) {
				dir_g = add_edges(dir_g, c(j, i))
			}
		}
		if (length(nn_out) > 0) {
			for (j in nn_out) {
				dir_g = add_edges(dir_g, c(i, j))
			}
		}
	}
	return(dir_g)
}

