#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>

using namespace Rcpp;

// ----------------------------
// Helper: Anderson-Darling A2
// ----------------------------
double compA2(const std::vector<double>& y) {
    int n = y.size();
    std::vector<double> z(y);
    double mean_z = std::accumulate(z.begin(), z.end(), 0.0)/n;
    double sd_z = 0.0;
    for(double v : z) sd_z += (v-mean_z)*(v-mean_z);
    sd_z = std::sqrt(sd_z/n);

    for(int i=0;i<n;i++) z[i] = (z[i]-mean_z)/sd_z;
    std::sort(z.begin(), z.end());

    const double eps = 1e-15;
    std::vector<double> Fi(n);
    for(int i=0;i<n;i++) {
        double p = 0.5*(1.0 + std::erf(z[i]/std::sqrt(2.0)));
        Fi[i] = std::min(std::max(p, eps), 1.0-eps);
    }

    double sum_term = 0.0;
    for(int i=0;i<n;i++) {
        sum_term += (2.0*(i+1)-1.0)*(std::log(Fi[i]) + std::log(1.0-Fi[n-1-i]));
    }

    double A2 = -n - sum_term/n;
    return A2 * (1.0 + 0.75/n + 2.25/(n*n));
}

// ----------------------------
// Helper: g-and-h transform
// ----------------------------
void gh_transform(const std::vector<double>& x, double g, double h, std::vector<double>& y) {
    int n = x.size();
    y.resize(n);
    for(int i=0;i<n;i++) {
        if(g != 0.0) y[i] = (std::exp(g*x[i])-1.0)/g * std::exp(h*x[i]*x[i]/2.0);
        else y[i] = x[i] * std::exp(h*x[i]*x[i]/2.0);
    }
}

// ----------------------------
// Objective function
// ----------------------------
double objfun(const std::vector<double>& x, double g, double h) {
    std::vector<double> y;
    gh_transform(x, g, h, y);
    return compA2(y);
}

// ----------------------------
// L-BFGS-B like optimizer (simplified)
// ----------------------------
void optimize2D(const std::vector<double>& x,
                double& g, double& h,
                double gmin, double gmax,
                double hmin, double hmax,
                int maxiter=100, double tol=1e-8) {

    double alpha = 1.0;    // line search step
    double beta = 0.5;     // backtracking factor
    double sigma = 1e-4;   // Armijo parameter

    for(int iter=0;iter<maxiter;iter++) {
        // Numerical gradients
        double eps = 1e-8;
        double f0 = objfun(x, g, h);

        double f_g = objfun(x, g+eps, h);
        double grad_g = (f_g - f0)/eps;

        double f_h = objfun(x, g, h+eps);
        double grad_h = (f_h - f0)/eps;

        // Check convergence
        if(std::abs(grad_g)<tol && std::abs(grad_h)<tol) break;

        // Gradient descent step with simple backtracking
        double step = alpha;
        double g_new = g - step*grad_g;
        double h_new = h - step*grad_h;

        // Project into bounds
        g_new = std::min(std::max(g_new, gmin), gmax);
        h_new = std::min(std::max(h_new, hmin), hmax);

        double f_new = objfun(x, g_new, h_new);

        // Armijo line search
        while(f_new > f0 - sigma*step*(grad_g*grad_g + grad_h*grad_h) && step>1e-12) {
            step *= beta;
            g_new = g - step*grad_g;
            h_new = h - step*grad_h;
            g_new = std::min(std::max(g_new, gmin), gmax);
            h_new = std::min(std::max(h_new, hmin), hmax);
            f_new = objfun(x, g_new, h_new);
        }

        g = g_new;
        h = h_new;
    }
}

// ----------------------------
// [[Rcpp::export]]
List osktfast_cpp(std::vector<double> x,
                  double g0=0.1, double h0=0.1,
                  double gmin=-1.0, double gmax=1.0,
                  double hmin=0.0, double hmax=0.5,
                  int maxiter=200) {

    // Standardize input
    double mean_x = std::accumulate(x.begin(), x.end(), 0.0)/x.size();
    double sd_x = 0.0;
    for(double v : x) sd_x += (v-mean_x)*(v-mean_x);
    sd_x = std::sqrt(sd_x/x.size());
    for(double& v : x) v = (v-mean_x)/sd_x;

    double g = g0, h = h0;
    optimize2D(x, g, h, gmin, gmax, hmin, hmax, maxiter);

    std::vector<double> y;
    gh_transform(x, g, h, y);

    return List::create(
        Named("transformed") = y,
        Named("g") = g,
        Named("h") = h,
        Named("value") = objfun(x, g, h)
    );
}
