# Back-transformation function for OSKT
backoskt <- function(Z, X_mean, X_sd, g, h, 
                     tol = 1e-10, maxiter = 1e6,
                     method = c("ur", "nr")) {
  
  method <- match.arg(method)
  
  t_start <- proc.time()[3]   # start timing
  
  # Newton-Raphson version (fast)
  # -------------------------------------------------
  if (method == "nr") {
      invert_single_nr <- function(z){
        xs <- z
        for(iter in 1:100){
          if(g != 0){
            Tx <- ((exp(g*xs)-1)/g)*exp(0.5*h*xs^2)
          } else {
            Tx <- xs*exp(0.5*h*xs^2)
          }
          fx <- Tx - z
          if(g != 0){
            dTx <- exp(0.5*h*xs^2)*exp(g*xs) + Tx*(h*xs)
          } else {
            dTx <- exp(0.5*h*xs^2)*(1 + h*xs^2)
          }
          step <- fx/dTx
          xs_new <- xs - step
          if(!is.finite(xs_new) || !is.finite(step)) {
            xs_new <- xs
            break
          }
          if(abs(xs_new - xs) < tol) return(xs_new)
          xs <- xs_new
        }
        return(xs)
     } 
     X_s <- sapply(Z, invert_single_nr)
  }
  
  # Uniroot version
  # -------------------------------------------------
  if (method == "ur") {
    
    invert_single_uniroot <- function(z) {
      f <- function(xs) {
        if (g != 0) {
          (((exp(g * xs) - 1) / g) * exp(0.5 * h * xs^2)) - z
        } else {
          xs * exp(0.5 * h * xs^2) - z
        }
      }
      
      lower <- -500
      upper <- 500
      
      root <- tryCatch(
        uniroot(f, lower = lower, upper = upper, tol = tol, maxiter = maxiter)$root,
        error = function(e) NA
      )
      return(root)
    }
    
    X_s <- sapply(Z, invert_single_uniroot)
  }
  
  X_orig <- X_s * X_sd + X_mean
  t_end <- proc.time()[3]
  
  return(list(
    X_orig = X_orig,
    X_s    = X_s,
    ct = t_end - t_start,
    method = method
  ))
}

