# Fast Inverse OSKT Transformation
backosktfast <- function(Z,
                         X_mean,
                         X_sd,
                         g,
                         h,
                         method = "auto",
                         tol = 1e-10,
                         maxiter_nr = 1000,
                         maxiter_brent = 2000) {

  ##---- Validate 'method' ----
  method <- match.arg(method, choices = c("auto", "nr", "brent"))

  ##---- Input validation ----
  if (!is.numeric(Z)) {
    stop("Z must be a numeric vector")
  }

  if (missing(X_mean) || !is.numeric(X_mean) || length(X_mean) != 1) {
    stop("X_mean must be provided as a single numeric value")
  }

  if (missing(X_sd) || !is.numeric(X_sd) || length(X_sd) != 1 || X_sd <= 0) {
    stop("X_sd must be provided as a single positive numeric value")
  }

  if (missing(g) || !is.numeric(g) || length(g) != 1) {
    stop("g must be provided as a single numeric value")
  }

  if (missing(h) || !is.numeric(h) || length(h) != 1) {
    stop("h must be provided as a single numeric value")
  }
  if (h < 0) {
    stop("h must be non-negative (h >= 0)")
  }

  if (!is.numeric(tol) || length(tol) != 1 || tol <= 0) {
    stop("tol must be a single positive numeric value")
  }

  if (!is.numeric(maxiter_nr) || length(maxiter_nr) != 1 ||
      maxiter_nr < 1 || maxiter_nr != as.integer(maxiter_nr)) {
    stop("maxiter_nr must be a positive integer")
  }

  if (!is.numeric(maxiter_brent) || length(maxiter_brent) != 1 ||
      maxiter_brent < 1 || maxiter_brent != as.integer(maxiter_brent)) {
    stop("maxiter_brent must be a positive integer")
  }


  ##---- Dispatch to C++ based on method ----
  if (method == "auto") {

    cpp_res <- backoskt_auto_cpp(
      Z = Z,
      X_mean = X_mean,
      X_sd = X_sd,
      g = g,
      h = h,
      tol = tol,
      maxiter_nr = maxiter_nr,
      maxiter_brent = maxiter_brent
    )

    # Convert numeric codes to meaningful character labels
    method_names <- character(length(cpp_res$method_used))
    method_names[cpp_res$method_used == 0] <- "failed"
    method_names[cpp_res$method_used == 1] <- "auto-nr"
    method_names[cpp_res$method_used == 2] <- "auto-brent"

    return(list(
      X_orig = cpp_res$X_orig,
      method_used = method_names
    ))

  } else if (method == "nr") {

    cpp_res <- backoskt_nr_cpp(
      Z = Z,
      X_mean = X_mean,
      X_sd = X_sd,
      g = g,
      h = h,
      tol = tol,
      maxiter = maxiter_nr
    )

    # Interpretation of method_used: 1=NR succeeded, 0=failed
    method_names <- ifelse(cpp_res$method_used == 1, "nr", "failed")

    return(list(
      X_orig = cpp_res$X_orig,
      method_used = method_names
    ))

  } else { # method == "brent"

    cpp_res <- backoskt_uniroot_cpp(
      Z = Z,
      X_mean = X_mean,
      X_sd = X_sd,
      g = g,
      h = h,
      tol = tol,
      maxiter = maxiter_brent
    )

    # Interpretation: 2=Brent succeeded, 0=failed
    method_names <- ifelse(cpp_res$method_used == 2, "brent", "failed")

    return(list(
      X_orig = cpp_res$X_orig,
      method_used = method_names
    ))
  }
}
