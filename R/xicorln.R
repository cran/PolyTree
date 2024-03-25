#' This function computes the xi correlation coefficient.
#' @param xvec The vector of x values.
#' @param yvec The vector of y values.

xicorln = function (xvec, yvec) {
    n <- length(xvec)
    PI <- rank(xvec, ties.method = "random")
    fr <- rank(yvec, ties.method = "max")/n
    gr <- rank((-yvec), ties.method = "max")/n
    ord <- order(PI)
    fr <- fr[ord]
    A1 <- sum(abs(fr[1:(n - 1)] - fr[2:n]))/(2 * n)
    CU <- mean(gr * (1 - gr))
    xi <- 1 - A1/CU
    return(xi)
}
