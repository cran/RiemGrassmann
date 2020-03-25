#' Pairwise Distance for Data on Grassmann Manifold
#' 
#' For data on grassmann manifold \eqn{x_1,x_2,\ldots,x_N \in Gr(k,n)}, compute pairwise distances \eqn{d(x_i,x_j)} via several metrics. The distance type \code{"intrinsic"} corresponds to 
#' geodesic distance while \code{"extrinsic"} is equivalent to \code{"chordal"} distance.
#' 
#' @param input either an array of size \eqn{(n\times k\times N)} or a list of length \eqn{N} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param type type of distance measure. Name of each type is \emph{Case Insensitive} and \emph{hyphen} can be omitted.
#' @param as.dist a logical; \code{TRUE} to return a \code{\link[stats]{dist}} object or \code{FALSE} to return an \eqn{(N\times N)} symmetric matrix.
#' @param useR a logical; \code{TRUE} to use R computations while \code{FALSE} goes everything in C++.
#' 
#' @return a \code{\link[stats]{dist}} object or \eqn{(N\times N)} symmetric matrix depending on \code{as.dist}.
#' 
#' @examples 
#' ## generate a dataset with two types of Grassmann elements
#' #  group1 : first four columns of (8x8) identity matrix + noise
#' #  group2 : last  four columns of (8x8) identity matrix + noise
#' 
#' mydata = list()
#' sdval  = 0.25
#' diag8  = diag(8)
#' for (i in 1:10){
#'   mydata[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' for (i in 11:20){
#'   mydata[[i]] = qr.Q(qr(diag8[,5:8] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' 
#' ## try 'intrinsic' distance using C++ implementation
#' dmat = gr.pdist(mydata, type="intrinsic", as.dist=FALSE)
#' opar = par(no.readonly=TRUE)
#' par(pty="s")
#' image(dmat, main="intrinsic distance")
#' par(opar)
#' 
#' \donttest{
#' ## compute and visualize distances for all types
#' #  we will iterate over all measures
#' alltypes = c("intrinsic","extrinsic","asimov","binet-cauchy",
#' "chordal","fubini-study","martin","procrustes","projection","spectral")
#' ntypes   = length(alltypes)
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (i in 1:ntypes){
#'   dmat = gr.pdist(mydata, type=alltypes[i], as.dist=FALSE)
#'   image(dmat[,20:1], main=alltypes[i])
#' }
#' par(opar)
#' }
#' 
#' @author Kisung You
#' @export
gr.pdist <- function(input, type=c("Intrinsic","Extrinsic","Asimov","Binet-Cauchy",
                               "Chordal","Fubini-Study","Martin","Procrustes",
                               "Projection","Spectral"), as.dist=TRUE, useR=FALSE){
  ############################################################
  # Preprocessing
  x        = RiemBase::riemfactory(return_gr(input), name="grassmann")
  alltypes = c("intrinsic","extrinsic","asimov","binetcauchy",
               "chordal","fubinistudy","martin","procrustes",
               "projection","spectral")
  if (missing(type)){
    mytype = "intrinsic"
  } else {
    mytype = match.arg(gsub("[-]","",tolower(type)),alltypes)
  }
  retdist = as.logical(as.dist)
  myuseR  = as.logical(useR)
  
  ############################################################
  # Computation
  if ((all(mytype=="intrinsic"))||(all(mytype=="extrinsic"))){
    output = RiemBaseExt::rstat.pdist(x, type=mytype, as.dist=FALSE)
  } else {
    if (myuseR){
      output = gr.pdist.others(x, method=mytype)  
    } else {
      output = cpp_pdist(x$data, mytype)
    }
  }
  
  ############################################################
  # Return
  if (retdist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}


# auxiliary here ----------------------------------------------------------
#' @keywords internal
#' @noRd
gr.pdist.others <- function(input, method){
  grdat = input$data
  N     = length(grdat)
  
  output = array(0,c(N,N))
  for (i in 1:(N-1)){
    Xi = as.matrix(grdat[[i]])
    for (j in (i+1):N){
      Xj = as.matrix(grdat[[j]])
      output[i,j] <- output[j,i] <- cdist_distance(Xi,Xj,method)
    }
  }
  return(output)
}
#' @keywords internal
#' @noRd
gr.pdist.nocheck <- function(input, mytype, asdist=TRUE){
  if ((all(mytype=="intrinsic"))||(all(mytype=="extrinsic"))){
    output = RiemBaseExt::rstat.pdist(input, type=mytype, as.dist=FALSE)
  } else {
    output = gr.pdist.others(input, method=mytype)
  }
  if (asdist){
    return(stats::as.dist(output))
  } else {
    return(output)
  }
}

# ## difference between C and R
# mydata = list()
# sdval  = 0.25
# diag8  = diag(8)
# for (i in 1:10){
#   mydata[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
# }
# for (i in 11:20){
#   mydata[[i]] = qr.Q(qr(diag8[,5:8] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
# }
# 
# #  we will iterate over all measures
# alltypes = c("intrinsic","extrinsic","asimov","binet-cauchy",
# "chordal","fubini-study","martin","procrustes","projection","spectral")
# ntypes   = length(alltypes)
# for (i in 1:ntypes){
#   dist.R = gr.pdist(mydata, type=alltypes[i], useR=TRUE, as.dist=FALSE)
#   dist.C = gr.pdist(mydata, type=alltypes[i], useR=FALSE, as.dist=FALSE)
#   print(paste0(alltypes[i]," difference=",norm(dist.R-dist.C,"F")))
# }
# 
# library(microbenchmark)
# microbenchmark::microbenchmark(
#   "R" = gr.pdist(mydata, type=alltypes[i], useR=TRUE, as.dist=FALSE),
#   "C" = gr.pdist(mydata, type=alltypes[i], useR=FALSE, as.dist=FALSE)
# )