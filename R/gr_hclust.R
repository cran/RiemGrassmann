#' Hierarchical Agglomerative Clustering on Grassmann Manifold
#' 
#' Given the \code{type} of distance measure and agglomeration scheme \code{method}, \code{gr.hclust} performs hierarchical clustering on 
#' Grassmann manifold using \pkg{fastcluster} package, which returns the same object as \pkg{stats} package's implementation while providing more efficient computation. 
#' See \code{\link[fastcluster]{hclust}} for more details.
#' 
#' @param input either an array of size \eqn{(n\times k\times N)} or a list of length \eqn{N} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param type type of distance measure. measure. Name of each type is \emph{Case Insensitive} and \emph{hyphen} can be omitted.
#' @param method he agglomeration method to be used. This must be (an unambiguous abbreviation of) one of \code{"single"}, \code{"complete"}, \code{"average"}, \code{"mcquitty"}, \code{"ward.D"}, \code{"ward.D2"}, \code{"centroid"} or \code{"median"}.
#' @param members \code{NULL} or a vector whose length equals the number of observations. See \code{\link[stats]{hclust}} for details.
#' 
#' @return an object of class \code{hclust}. See \code{\link[stats]{hclust}} for details. 
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
#' ## try hierarchical clustering with "intrinsic" distance
#' opar <- par(no.readonly=TRUE)
#' hint <- gr.hclust(mydata, type="intrinsic")
#' plot(hint, main="intrinsic+single")
#' par(opar)
#' 
#' \donttest{
#' ## do hierarchical clustering with different distance measures
#' alltypes = c("intrinsic","extrinsic","asimov","binet-cauchy",
#' "chordal","fubini-study","martin","procrustes","projection","spectral")
#' ntypes   = length(alltypes)
#' 
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (i in 1:ntypes){
#'   hout = gr.hclust(mydata, type=alltypes[i])
#'   plot(hout, main=paste0("hclust::",alltypes[i]))
#' }
#' par(opar)
#' }
#' 
#' @author Kisung You
#' @export
gr.hclust <- function(input, type=c("Intrinsic","Extrinsic","Asimov","Binet-Cauchy",
                                    "Chordal","Fubini-Study","Martin","Procrustes",
                                    "Projection","Spectral"),
                      method = c("single", "complete", "average", "mcquitty", "ward.D", "ward.D2",
                                 "centroid", "median"),
                      members=NULL){
  ############################################################
  # Preprocessing
  x        = RiemBase::riemfactory(return_gr(input), name="grassmann")
  alltypes = c("intrinsic","extrinsic","asimov","binetcauchy",
               "chordal","fubinistudy","martin","procrustes",
               "projection","spectral")
  if (missing(type)){
    mydtype = "intrinsic"
  } else {
    mydtype = match.arg(gsub("[-]","",tolower(type)),alltypes)
  }
  
  ############################################################
  # Compute Distance and Apply Hclust
  pdmat = gr.pdist.nocheck(x, mydtype, asdist=TRUE)
  hcout = RiemBaseExt::rclust.hclust(pdmat)
  return(hcout)
}