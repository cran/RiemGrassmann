#' Pairwise Distance for Two Sets Data on Grassmann Manifold
#' 
#' For data on grassmann manifold \eqn{x_1,x_2,\ldots,x_M \in Gr(k,n)} and 
#' \eqn{y_1,y_2,\ldots,y_N \in Gr(k,n)}, compute pairwise distances \eqn{d(x_i,y_j)} via several metrics. The distance type \code{"intrinsic"} corresponds to 
#' geodesic distance while \code{"extrinsic"} is equivalent to \code{"chordal"} distance.
#' 
#' 
#' @param input1 either an array of size \eqn{(n\times k\times M)} or a list of length \eqn{M} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param input2 either an array of size \eqn{(n\times k\times N)} or a list of length \eqn{N} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param type type of distance measure. Name of each type is \emph{Case Insensitive} and \emph{hyphen} can be omitted.
#' @param useR a logical; \code{TRUE} to use R computations while \code{FALSE} goes everything in C++.
#' 
#' @return an \eqn{(M\times N)} matrix of pairwise distances.
#' 
#' @examples 
#' ## generate a dataset with two types of Grassmann elements
#' #  group1 : first four columns of (8x8) identity matrix + noise
#' #  group2 : last  four columns of (8x8) identity matrix + noise
#' 
#' mydata1 = list()
#' mydata2 = list()
#' sdval   = 0.25
#' diag8   = diag(8)
#' for (i in 1:10){
#'   mydata1[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' for (i in 1:10){
#'   mydata2[[i]] = qr.Q(qr(diag8[,5:8] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
#' }
#' 
#' 
#' \donttest{
#' ## compute and visualize distances for all types
#' #  we will iterate over all measures
#' alltypes = c("intrinsic","extrinsic","asimov","binetcauchy",
#' "chordal","fubinistudy","martin","procrustes","projection","spectral")
#' ntypes   = length(alltypes)
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (i in 1:ntypes){
#'   dmat = gr.pdist2(mydata1, mydata2, type=alltypes[i])
#'   image(dmat, main=alltypes[i])
#' }
#' par(opar)
#' }
#' 
#' @author Kisung You
#' @export
gr.pdist2 <- function(input1, input2, type=c("Intrinsic","Extrinsic","Asimov","Binet-Cauchy",
                                             "Chordal","Fubini-Study","Martin","Procrustes",
                                             "Projection","Spectral"), useR=FALSE){
  ############################################################
  # Preprocessing
  input1 = RiemBase::riemfactory(return_gr(input1), name="grassmann")
  input2 = RiemBase::riemfactory(return_gr(input2), name="grassmann")
  alltypes = c("intrinsic","extrinsic","asimov","binetcauchy",
               "chordal","fubinistudy","martin","procrustes",
               "projection","spectral")
  if (missing(type)){
    mytype = "intrinsic"
  } else {
    mytype = match.arg(gsub("[-]","",tolower(type)),alltypes)
  }
  myuseR = as.logical(useR)

  ############################################################
  # Computation
  if ((all(mytype=="intrinsic"))||(all(mytype=="extrinsic"))){
    output = RiemBaseExt::rstat.pdist2(input1, input2, type=mytype)
  } else {
    if (myuseR){
      output = gr.pdist2.others(input1, input2, method=mytype)  
    } else {
      output = cpp_pdist2(input1$data, input2$data, mytype)
    }
  }
  
  ############################################################
  # Return
  return(output)
}

# auxiliary here ----------------------------------------------------------
#' @keywords internal
#' @noRd
gr.pdist2.others <- function(input1, input2, method){
  gdat1 = input1$data
  gdat2 = input2$data
  
  N = length(gdat1)
  M = length(gdat2)
  
  output = array(0,c(N,M))
  for (n in 1:N){
    Xn = as.matrix(gdat1[[n]])
    for (m in 1:M){
      Xm = as.matrix(gdat2[[m]])
      output[n,m] <- cdist_distance(Xn,Xm,method)
    }
  }
  return(output)
}

# ## personal test
# mydata1 = list()
# mydata2 = list()
# sdval   = 0.25
# diag8   = diag(8)
# for (i in 1:10){
#   mydata1[[i]] = qr.Q(qr(diag8[,1:4] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
# }
# for (i in 1:10){
#   mydata2[[i]] = qr.Q(qr(diag8[,5:8] + matrix(rnorm(8*4,sd=sdval),ncol=4)))
# }
# 
# alltypes = c("intrinsic","extrinsic","asimov","binet-cauchy",
# "chordal","fubini-study","martin","procrustes","projection","spectral")
# ntypes   = length(alltypes)
# for (i in 1:ntypes){
#   dist.R = gr.pdist2(mydata1, mydata2, type=alltypes[i], useR=TRUE)
#   dist.C = gr.pdist2(mydata1, mydata2, type=alltypes[i], useR=FALSE)
#   print(paste0(alltypes[i]," difference=",norm(dist.R-dist.C,"F")))
# }
# 
# library(microbenchmark)
# microbenchmark::microbenchmark(
#   "R" = gr.pdist2(mydata1, mydata2, type=alltypes[i], useR=TRUE),
#   "C" = gr.pdist2(mydata1, mydata2, type=alltypes[i], useR=FALSE)
# )