#include <Rcpp.h>
#include <cmath>
#include <limits>

using namespace Rcpp;

// -----------------------------------------------------------
// OSKT forward function and derivative
// -----------------------------------------------------------
inline double oskt_fx(double xs, double z, double g, double h) {
  double exp_h = std::exp(0.5 * h * xs * xs);
  if (std::abs(g) > 1e-8) {
    return ((std::exp(g * xs) - 1.0) / g) * exp_h - z;
  } else {
    return xs * exp_h - z;
  }
}

inline double oskt_dfx(double xs, double g, double h) {
  double exp_h = std::exp(0.5 * h * xs * xs);
  if (std::abs(g) > 1e-8) {
    double Tx = ((std::exp(g * xs) - 1.0) / g) * exp_h;
    return exp_h * std::exp(g * xs) + Tx * (h * xs);
  } else {
    return exp_h * (1.0 + h * xs * xs);
  }
}

// -----------------------------------------------------------
// 1) Newton–Raphson (fast, derivative-based root-finding)
// -----------------------------------------------------------
// [[Rcpp::export]]
List backoskt_nr_cpp(
    NumericVector Z,
    double X_mean,
    double X_sd,
    double g,
    double h,
    double tol = 1e-10,
    int maxiter = 1000
) {
  // Parameter validation
  if (X_sd <= 0.0) {
    stop("X_sd must be positive (X_sd > 0)");
  }
  if (h < 0.0) {
    stop("Parameter 'h' must be non-negative (h >= 0)");
  }
  
  int n = Z.size();
  NumericVector Xs(n);
  IntegerVector method_used(n);  // 0=failed, 1=nr

  for (int i = 0; i < n; i++) {

    double z = Z[i];
    if (!R_finite(z)) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
      continue;
    }

    // Robust initial guess
    double xs = (std::abs(z) < 1.0)
                  ? z
                  : std::copysign(std::log1p(std::abs(z)), z);

    bool converged = false;

    for (int iter = 0; iter < maxiter; iter++) {
      double fx  = oskt_fx(xs, z, g, h);
      double dfx = oskt_dfx(xs, g, h);

      if (!R_finite(fx) || !R_finite(dfx) || std::abs(dfx) < 1e-12)
        break;

      double step = fx / dfx;
      double xs_new = xs - step;

      if (!R_finite(xs_new))
        break;

      // Dual convergence criterion with residual check
      if (std::abs(xs_new - xs) < tol && std::abs(fx) < tol) {
        xs = xs_new;
        // Extra validation: verify final residual
        double final_fx = oskt_fx(xs, z, g, h);
        if (std::abs(final_fx) < 10.0 * tol) {
          converged = true;
          break;
        }
      }

      xs = xs_new;
    }

    if (converged) {
      Xs[i] = xs * X_sd + X_mean;
      method_used[i] = 1;
    } else {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
    }
  }

  return List::create(
    Named("X_orig") = Xs,
    Named("method_used") = method_used
  );
}

// -----------------------------------------------------------
// 2) Brent–Dekker (robust, bracketing root-finding)
// -----------------------------------------------------------
// [[Rcpp::export]]
List backoskt_uniroot_cpp(
    NumericVector Z,
    double X_mean,
    double X_sd,
    double g,
    double h,
    double tol = 1e-10,
    int maxiter = 2000
) {
  // Parameter validation
  if (X_sd <= 0.0) {
    stop("X_sd must be positive (X_sd > 0)");
  }
  if (h < 0.0) {
    stop("Parameter 'h' must be non-negative (h >= 0)");
  }
  
  int n = Z.size();
  NumericVector Xs(n);
  IntegerVector method_used(n);  // 0=failed, 2=brent

  for (int i = 0; i < n; i++) {

    double z = Z[i];
    if (!R_finite(z)) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
      continue;
    }

    // Adaptive initial bracket based on z
    double scale = std::max(10.0, 
                           std::sqrt(2.0 / std::max(h, 1e-8)) * 
                           std::log1p(std::abs(z)));
    double a = -scale, b = scale;
    double fa = oskt_fx(a, z, g, h);
    double fb = oskt_fx(b, z, g, h);

    // Expand bracket if needed
    int expand_attempts = 0;
    while (fa * fb > 0 && expand_attempts < 20) {
      if (std::abs(fa) < std::abs(fb)) {
        a = a - (b - a);
        fa = oskt_fx(a, z, g, h);
      } else {
        b = b + (b - a);
        fb = oskt_fx(b, z, g, h);
      }
      expand_attempts++;
    }

    if (fa * fb > 0) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
      continue;
    }

    // Ensure sign convention: fa < 0, fb > 0
    if (fa > 0) {
      std::swap(a, b);
      std::swap(fa, fb);
    }

    double c = a, fc = fa;
    double d = 0.0;
    bool mflag = true;

    for (int iter = 0; iter < maxiter; iter++) {

      // Convergence check
      if (std::abs(b - a) < tol || std::abs(fb) < tol) {
        break;
      }

      double s;

      // Inverse quadratic interpolation or secant
      if (fa != fc && fb != fc) {
        // Inverse quadratic interpolation
        s = (a * fb * fc) / ((fa - fb) * (fa - fc))
          + (b * fa * fc) / ((fb - fa) * (fb - fc))
          + (c * fa * fb) / ((fc - fa) * (fc - fb));
      } else {
        // Secant method
        s = b - fb * (b - a) / (fb - fa);
      }

      // Bisection safety conditions
      bool condition1 = (s < (3.0 * a + b) / 4.0) || (s > b);
      bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0);
      bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0);
      bool condition4 = mflag && (std::abs(b - c) < tol);
      bool condition5 = !mflag && (std::abs(c - d) < tol);

      if (condition1 || condition2 || condition3 || condition4 || condition5) {
        s = (a + b) / 2.0;  // Use bisection
        mflag = true;
      } else {
        mflag = false;
      }

      double fs = oskt_fx(s, z, g, h);
      d = c;
      c = b;
      fc = fb;

      // Update brackets
      if (fa * fs < 0) {
        b = s;
        fb = fs;
      } else {
        a = s;
        fa = fs;
      }

      // Keep best estimate in 'b'
      if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
      }
    }

    // Final residual check
    double final_residual = std::abs(oskt_fx(b, z, g, h));
    if (final_residual > 10.0 * tol) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
    } else {
      Xs[i] = b * X_sd + X_mean;
      method_used[i] = 2;
    }
  }

  return List::create(
    Named("X_orig") = Xs,
    Named("method_used") = method_used
  );
}

// -----------------------------------------------------------
// 3) AUTO: NR with Brent fallback (recommended)
// -----------------------------------------------------------
// [[Rcpp::export]]
List backoskt_auto_cpp(
    NumericVector Z,
    double X_mean,
    double X_sd,
    double g,
    double h,
    double tol = 1e-10,
    int maxiter_nr = 1000,
    int maxiter_brent = 2000
) {
  // Parameter validation
  if (X_sd <= 0.0) {
    stop("X_sd must be positive (X_sd > 0)");
  }
  if (h < 0.0) {
    stop("Parameter 'h' must be non-negative (h >= 0)");
  }
  
  int n = Z.size();
  NumericVector Xs(n);
  IntegerVector method_used(n);  // 0=failed, 1=NR, 2=Brent

  for (int i = 0; i < n; i++) {

    double z = Z[i];
    if (!R_finite(z)) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
      continue;
    }

    // ========================================================
    // PHASE 1: Newton-Raphson
    // ========================================================
    double xs = (std::abs(z) < 1.0)
                  ? z
                  : std::copysign(std::log1p(std::abs(z)), z);
    
    bool nr_converged = false;
    
    for (int iter = 0; iter < maxiter_nr; iter++) {
      double fx  = oskt_fx(xs, z, g, h);
      double dfx = oskt_dfx(xs, g, h);

      if (!R_finite(fx) || !R_finite(dfx) || std::abs(dfx) < 1e-12)
        break;

      double step = fx / dfx;
      double xs_new = xs - step;

      if (!R_finite(xs_new))
        break;

      if (std::abs(xs_new - xs) < tol && std::abs(fx) < tol) {
        xs = xs_new;
        // Extra validation: verify final residual
        double final_fx = oskt_fx(xs, z, g, h);
        if (std::abs(final_fx) < 10.0 * tol) {
          nr_converged = true;
        }
        break;
      }

      xs = xs_new;
    }

    if (nr_converged) {
      Xs[i] = xs * X_sd + X_mean;
      method_used[i] = 1;
      continue;
    }

    // ========================================================
    // PHASE 2: Brent-Dekker fallback
    // ========================================================
    double scale = std::max(10.0, 
                           std::sqrt(2.0 / std::max(h, 1e-8)) * 
                           std::log1p(std::abs(z)));
    double a = -scale, b = scale;
    double fa = oskt_fx(a, z, g, h);
    double fb = oskt_fx(b, z, g, h);

    int expand_attempts = 0;
    while (fa * fb > 0 && expand_attempts < 20) {
      if (std::abs(fa) < std::abs(fb)) {
        a = a - (b - a);
        fa = oskt_fx(a, z, g, h);
      } else {
        b = b + (b - a);
        fb = oskt_fx(b, z, g, h);
      }
      expand_attempts++;
    }

    if (fa * fb > 0) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
      continue;
    }

    if (fa > 0) {
      std::swap(a, b);
      std::swap(fa, fb);
    }

    double c = a, fc = fa;
    double d = 0.0;
    bool mflag = true;

    for (int iter = 0; iter < maxiter_brent; iter++) {
      if (std::abs(b - a) < tol || std::abs(fb) < tol) {
        break;
      }

      double s;

      if (fa != fc && fb != fc) {
        s = (a * fb * fc) / ((fa - fb) * (fa - fc))
          + (b * fa * fc) / ((fb - fa) * (fb - fc))
          + (c * fa * fb) / ((fc - fa) * (fc - fb));
      } else {
        s = b - fb * (b - a) / (fb - fa);
      }

      bool condition1 = (s < (3.0 * a + b) / 4.0) || (s > b);
      bool condition2 = mflag && (std::abs(s - b) >= std::abs(b - c) / 2.0);
      bool condition3 = !mflag && (std::abs(s - b) >= std::abs(c - d) / 2.0);
      bool condition4 = mflag && (std::abs(b - c) < tol);
      bool condition5 = !mflag && (std::abs(c - d) < tol);

      if (condition1 || condition2 || condition3 || condition4 || condition5) {
        s = (a + b) / 2.0;
        mflag = true;
      } else {
        mflag = false;
      }

      double fs = oskt_fx(s, z, g, h);
      d = c;
      c = b;
      fc = fb;

      if (fa * fs < 0) {
        b = s;
        fb = fs;
      } else {
        a = s;
        fa = fs;
      }

      if (std::abs(fa) < std::abs(fb)) {
        std::swap(a, b);
        std::swap(fa, fb);
      }
    }

    double final_residual = std::abs(oskt_fx(b, z, g, h));
    if (final_residual > 10.0 * tol) {
      Xs[i] = NA_REAL;
      method_used[i] = 0;
    } else {
      Xs[i] = b * X_sd + X_mean;
      method_used[i] = 2;
    }
  }

  return List::create(
    Named("X_orig") = Xs,
    Named("method_used") = method_used
  );
}