#' Inference, Learning, and Optimization on Grassmann Manifold
#' 
#' Grassmann manifold \eqn{Gr(k,n)} is the set of k-planes, or \eqn{k}-dimensional subspaces in \eqn{R^n}, 
#' which is indeed a compact Riemannian manifold. In this package, we use a convention that each element 
#' in \eqn{Gr(k,n)} is represented as an orthonormal basis (ONB) \eqn{X \in \mathbf{R}^{n\times k}} where
#' \deqn{X^\top X = I_k}. 
#' 
#' @docType package
#' @name package-RiemGrassmann
#' @importFrom RiemBase riemfactory
#' @importFrom RiemBaseExt rstat.frechet rstat.pdist rstat.pdist2 rclust.hclust rclust.kmedoids
#' @importFrom stats rnorm runif cov cmdscale
#' @importFrom utils packageVersion
#' @importFrom Rcpp evalCpp
#' @useDynLib RiemGrassmann
NULL
# https://en.wikipedia.org/wiki/List_of_probability_distributions
# pack <- "T4mle"
# path <- find.package(pack)
# system(paste(shQuote(file.path(R.home("bin"), "R")),
#              "CMD", "Rd2pdf", shQuote(path)))