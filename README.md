# OSKTNORM: Normalization via Optimized Skewness and Kurtosis of Non-Normal Data

## Introduction

The `osktnorm` package provides functions for normalizing non-normal data by optimizing skewness and kurtosis simultaneously. It can handle both right- and left-skewed distributions, 
and it includes functions for forward transformation (`osktfast`) and back-transformation (`backosktfast`). The package also provides normality test statistics to evaluate transformation performance.

---

## Install osktnorm

You can install the package from CRAN using:

```r
install.packages("osktnorm", repos="https://cloud.r-project.org", dep=TRUE)
```

After installation, load the package:

```r
library(osktnorm)
```

---

## Quick OSKT Normalization

### Generate example data

```r
set.seed(12)
x_orig <- rlnorm(300, mean=0, sd=0.5)  # Right-skewed data
```

### Apply OSKT

```r
res_oskt <- osktfast(x_orig)           # Forward transformation
x_trans <- res_oskt$transformed        # Transformed values
g_star <- res_oskt$g                    # Optimized skewness parameter
h_star <- res_oskt$h                    # Optimized kurtosis parameter
A2 <- res_oskt$value                    # Anderson-Darling statistic

head(x_trans, 5)
cat("Optimized skewness:", g_star, "Optimized kurtosis:", h_star, "A2:", A2, "\n")
```

### Visualization

```r
breaks <- pretty(range(c(x_orig, x_trans)), n = 25)
hist(x_orig, breaks=breaks, freq=FALSE, col=rgb(0.2,0.4,0.8,0.4), main="Before and After OSKT", xlab="Value")
lines(density(x_orig), col="blue")
hist(x_trans, breaks=breaks, freq=FALSE, col=rgb(0.8,0.3,0.3,0.4), add=TRUE)
lines(density(x_trans), col="red")
curve(dnorm(x), add=TRUE, lty=2, col="black")
legend("topleft", legend=c("Original","Transformed","Original Density","OSKT Density","Standard Normal"),
       col=c(rgb(0.2,0.4,0.8,0.6), rgb(0.8,0.3,0.3,0.6), "blue","red","black"),
       lty=c(1,1,1,1,2), lwd=c(10,10,2,2,2), bty="n")
```

---

## Back-transformation

Recover original values using `backosktfast`:

```r
res_back <- backosktfast(Z=x_trans, X_mean=mean(x_orig), X_sd=sd(x_orig),
                          g=g_star, h=h_star, method="auto")
x_recovered <- res_back$X_orig
all.equal(x_orig, x_recovered, tolerance=1e-6)
```

### Plot recovered vs original

```r
plot(x_orig, x_recovered, pch=16, col=rgb(0.2,0.2,0.7,0.4), xlab="Original", ylab="Recovered")
abline(0,1, col="red", lty=2)
```

---

## Normality Comparison with Other Methods

```r
# Load example left-skewed data
x_orig <- groupcompare::ghdist(n=300, A=0, B=1, g=-0.49, h=0)

# Apply transformations
x_bc   <- boxcox(x_orig, makepositive=TRUE)$transformed
x_yj   <- yeojohnson(x_orig)$transformed
x_oskt <- osktfast(x_orig)$transformed

# Normality statistics
get_stats <- function(x) {
  x <- x[is.finite(x)]
  c(Skew=mean((x-mean(x))^3)/sd(x)^3,
    Kurt=mean((x-mean(x))^4)/sd(x)^4-3,
    SW=shapiro.test(x)$p.value,
    ZA=zatest(x, nsim=100)$p.value,
    CVM=cvmtest(x)$p.value,
    PPM=unname(pearsonp(x)$statistic))
}

pval_table <- rbind(ORG=get_stats(x_orig),
                     BC=get_stats(x_bc),
                     YJ=get_stats(x_yj),
                     OSKT=get_stats(x_oskt))
round(pval_table,4)
```

---

## Citation

To cite the `osktnorm` package in publications, run:

```r
citation("osktnorm")
```
