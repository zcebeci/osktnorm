# High-speed OSKT Transformation with Anderson-Darling Minimization
osktfast <- function(x,
                     init_params = c(0.1, 0.1),
                     lower_bounds = c(-1, 0),
                     upper_bounds = c(1, 0.5),
                     maxiter = 200) {
  
  if(length(init_params) != 2) stop("init_params must be length 2 (g, h).")
  if(length(lower_bounds) != 2 || length(upper_bounds) != 2) stop("Bounds must be length 2.")
  
  g0 <- init_params[1]
  h0 <- init_params[2]
  
  gmin <- lower_bounds[1]
  hmin <- lower_bounds[2]
  
  gmax <- upper_bounds[1]
  hmax <- upper_bounds[2]
  
  # Call the C++ function
  res <- osktfast_cpp(x,
                      g0 = g0, h0 = h0,
                      gmin = gmin, gmax = gmax,
                      hmin = hmin, hmax = hmax,
                      maxiter = maxiter)
  
  return(res)
}
