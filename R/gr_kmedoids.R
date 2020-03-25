#' k-Medoids Clustering on Grassmann Manifold
#' 
#' k-Medoids algorithm depends solely on the availability of concept that gives dissimilarity. We adopt \code{pam} algorithm 
#' from \pkg{cluster} package. See \code{\link[cluster]{pam}} for more details.
#' 
#' @param input either an array of size \eqn{(n\times k\times N)} or a list of length \eqn{N} whose elements are \eqn{(n\times k)} orthonormal basis (ONB) on Grassmann manifold.
#' @param k the number of clusters
#' @param type type of distance measure. measure. Name of each type is \emph{Case Insensitive} and \emph{hyphen} can be omitted.
#' 
#' @return an object of class \code{pam}. See \code{\link[cluster]{pam}} for details.
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
#' ## do k-medoids clustering with 'intrinsic' distance
#' #  First, apply MDS for visualization
#' dmat = gr.pdist(mydata, type="intrinsic")
#' embd = stats::cmdscale(dmat, k=2)
#' 
#' #  Run 'gr.kmedoids' with different numbers of clusters
#' grint2 = gr.kmedoids(mydata, type="intrinsic", k=2)$clustering
#' grint3 = gr.kmedoids(mydata, type="intrinsic", k=3)$clustering
#' grint4 = gr.kmedoids(mydata, type="intrinsic", k=4)$clustering
#' 
#' #  Let's visualize
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(1,3), pty="s")
#' plot(embd, pch=19, col=grint2, main="k=2")
#' plot(embd, pch=19, col=grint3, main="k=3")
#' plot(embd, pch=19, col=grint4, main="k=4")
#' par(opar)
#' 
#' \donttest{
#' ## perform k-medoids clustering with different distance measures
#' #  iterate over all distance measures
#' alltypes = c("intrinsic","extrinsic","asimov","binet-cauchy",
#' "chordal","fubini-study","martin","procrustes","projection","spectral")
#' ntypes   = length(alltypes)
#' labels   = list()
#' for (i in 1:ntypes){
#'   labels[[i]] = gr.kmedoids(mydata, k=2, type=alltypes[i])$clustering
#' }
#' 
#' ## visualize
#' #  1. find MDS scaling for each distance measure as well
#' embeds = list()
#' for (i in 1:ntypes){
#'   pdmat       = gr.pdist(mydata, type=alltypes[i])
#'   embeds[[i]] = stats::cmdscale(pdmat, k=2)
#' }
#' 
#' #  2. plot the clustering results
#' opar <- par(no.readonly=TRUE)
#' par(mfrow=c(2,5), pty="s")
#' for (i in 1:ntypes){
#'   pm = paste0("k-medoids::",alltypes[i])
#'   plot(embeds[[i]], col=labels[[i]], main=pm, pch=19)
#' }
#' par(opar)
#' }
#' 
#' @author Kisung You
#' @export
gr.kmedoids <- function(input, k=2,
                        type=c("Intrinsic","Extrinsic","Asimov","Binet-Cauchy",
                               "Chordal","Fubini-Study","Martin","Procrustes",
                               "Projection","Spectral")){
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
  myk = round(k)
  
  ############################################################
  # Compute Pairwise Distance and Run 'rclust.kmedoids'
  pdmat  = gr.pdist.nocheck(x, mydtype, asdist=TRUE)
  output = RiemBaseExt::rclust.kmedoids(pdmat, k=myk)
  return(output)
}
