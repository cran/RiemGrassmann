# computations
# (1) cdist_prangles : compute principal angles between two subspaces.
# (2) cdist_distance : compute distances other than intrinsic/extrinsic.

# (1) cdist_prangles
#' @keywords internal
cdist_prangles <- function(U,V){
  return(as.vector(base::acos((base::svd(t(U)%*%V)$d))))
}
# (2) cdist_distance
#' @keywords internal
cdist_distance <- function(U,V,method){
  # compute principal angles
  k = ncol(U)
  thetas = cdist_prangles(U,V) # in a descending order
  
  # switching cases via 'if' clauses
  if (all(method=="asimov")){
    return(thetas[k])
  } else if (all(method=="binetcauchy")){
    return(sqrt(1-(base::prod((cos(thetas))^2))))
  } else if (all(method=="chordal")){
    return(sqrt(sum(base::sin(thetas)^2)))
  } else if (all(method=="fubinistudy")){
    return(base::acos(base::prod(base::cos(thetas))))
  } else if (all(method=="martin")){
    return(sqrt(sum(-2*log(cos(thetas)))))
  } else if (all(method=="procrustes")){
    return(2*sqrt(sum((sin(thetas/2)^2))))
  } else if (all(method=="projection")){
    return(sin(thetas[k]))
  } else if (all(method=="spectral")){
    return(sin(thetas[k]/2)*2)
  }
}