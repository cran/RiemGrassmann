#' Fréchet Mean on Grassmann Manifold
#' 
#' For manifold-valued data, Fréchet mean is the solution of following cost function,
#' \deqn{\textrm{min}_x \sum_{i=1}^n \rho^2 (x, x_i),\quad x\in\mathcal{M}}
#' for a given data \eqn{\{x_i\}_{i=1}^n} and \eqn{\rho(x,y)} is the geodesic distance 
#' between two points on manifold \eqn{\mathcal{M}}. It uses a gradient descent method 
#' with a backtracking search rule for updating. 
#' 
#' @param x either an array of size \eqn{(n\times k\times N)} or a list of length \eqn{N} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param type type of geometry, either \code{"intrinsic"} or \code{"extrinsic"}.
#' @param eps stopping criterion for the norm of gradient.
#' @param parallel a flag for enabling parallel computation with OpenMP.
#' 
#' @return a named list containing
#' \describe{
#' \item{mu}{an estimated mean matrix for ONB of size \eqn{(n\times k)}.}
#' \item{variation}{Fréchet variation with the estimated mean.}
#' }
#' 
#' @examples 
#' ## generate a dataset with two types of Grassmann elements
#' #  first four columns of (8x8) identity matrix + noise
#' mydata = list()
#' sdval  = 0.1
#' diag8  = diag(8)
#' for (i in 1:10){
#'   mydata[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' 
#' ## compute two types of means
#' mean.int = gr.mean(mydata, type="intrinsic")
#' mean.ext = gr.mean(mydata, type="extrinsic")
#' 
#' ## visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,2))
#' image(mean.int$mu, main="intrinsic mean")
#' image(mean.ext$mu, main="extrinsic mean")
#' par(opar)
#' 
#' @author Kisung You
#' @export
gr.mean <- function(x, type=c("intrinsic","extrinsic"), eps=1e-6, parallel=FALSE){
  ############################################################
  # Preprocessing
  x      = RiemBase::riemfactory(return_gr(x), name="grassmann")
  mytype = match.arg(type)
  myeps      = as.double(eps)
  myparallel = as.logical(parallel)
  
  ############################################################
  # Computation
  output = RiemBaseExt::rstat.frechet(x, type=mytype, int.eps=myeps, parallel=myparallel)
  return(output)
}