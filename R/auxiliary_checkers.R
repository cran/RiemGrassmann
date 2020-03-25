## checkers
## (1) return_gr : return grassmann manifold as 3d object

# (1) return_gr
#' @keywords internal
return_gr <- function(x){
  wow = RiemBase::riemfactory(x, name="grassmann")
  N = length(wow$data)
  n = wow$size[1]
  p = wow$size[2]
  
  output = array(0,c(n,p,N))
  for (i in 1:N){
    output[,,i] = wow$data[[i]]
  }
  return(output)
}