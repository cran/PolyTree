#' This is the function that computes the skeletion tree from data. 
#' The input is a matrix x whose rows are the data vectors. 
#' The sample size n is the number of rows.
#' The number of variables p is the number of columns
#' The function outputs the skeleton tree g.
#' @param x The input data matrix.

condeptree = function(x) {
	p = ncol(x)
	n = nrow(x)
 	d = matrix(0, nrow = p, ncol = p)
	for (i in 1:(p-1)) {
		for (j in (i+1):p) {
			d[i,j] = xicorln(x[,i], x[,j])
			d[j,i] = xicorln(x[,j], x[,i])
		}
	}
	m = d
	maxval = max(d)
	for (i in 1:(p-1)) {
		for (j in (i+1):p) {
			u = d[,i]
			v = d[,j]
			w1 = setdiff(which(u >= d[j,i]), c(i,j))
			w2 = setdiff(which(v >= d[i,j]), c(i,j))
			if (length(intersect(w1,w2)) > 0) {
				m[i,j] = 0
				m[j,i] = 0
			}
			if (m[i,j] != 0 & m[j,i] != 0) {
				u1 = maxval + 1 - m[i,j]
				u2 = maxval + 1 - m[j,i]
				u = max(u1,u2)
				m[i,j] = u
				m[j,i] = u
			} 
		}
	}
	m = graph_from_adjacency_matrix(m, mode = "undirected", weighted = TRUE)
	g = mst(m)
	return(g)
}
