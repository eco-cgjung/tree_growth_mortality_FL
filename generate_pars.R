generate.pars <- function(p_op,pmin,pmax,d) {
  while(TRUE){
    rand <- runif(length(pmin))
    pnew <- p_op+(rand-0.5)*(pmax-pmin)/d
    if  (Reduce("&", pnew>pmin&pnew<pmax))
      break
  }
  pnew
}

generate.pars.cov <- function(p_op,pmin,pmax, covars) {
  while(TRUE){
    pnew <- rmvn(1, mu = p_op, sigma = covars,ncores = 12)
    if  (Reduce("&", pnew>pmin&pnew<pmax))
      break
  }
  pnew
}