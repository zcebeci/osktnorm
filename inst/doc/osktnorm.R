knitr::opts_chunk$set(
  collapse = TRUE,
  comment = " ",
  options("width" = 120),
  fig.width=7, fig.height=5
)

# install.packages("osktnorm", repos="https://cloud.r-project.org", dep=TRUE)

library(osktnorm)

set.seed(12)
x_orig <- rlnorm(300, mean=0, sd=0.5) # Generate right-skewed data

# Apply OSKT normality
res_oskt <- osktfast(x_orig) 

x_transformed <- res_oskt$transformed
head(x_transformed, 5)

g_star <- res_oskt$g
h_star <- res_oskt$h
A2 <- res_oskt$value

cat("Optimized skewness: ", g_star, "\n")
cat("Optimized kurtosis: ", h_star, "\n")
cat("Anderson-Darling statistic at the optimum: ", A2, "\n")

breaks <- pretty(range(c(x_orig, x_transformed)), n = 25)
h_orig  <- hist(x_orig, breaks = breaks, plot = FALSE)
h_trans <- hist(x_transformed, breaks = breaks, plot = FALSE)

d_orig  <- density(x_orig); d_trans <- density(x_transformed)

ymax <- max(c(h_orig$density, h_trans$density, d_orig$y,d_trans$y, dnorm(0)))
hist(x_orig, breaks = breaks, freq = FALSE, ylim = c(0, ymax * 1.05), 
     col = rgb(0.2, 0.4, 0.8, 0.4), border = "white", 
     main = "Before and After OSKT Transformation", xlab = "Value")
lines(d_orig, col = "blue", lwd = 2)

hist(x_transformed, breaks = breaks, freq = FALSE,
     col = rgb(0.8, 0.3, 0.3, 0.4), border = "white", add = TRUE)
lines(d_trans, col = "red", lwd = 2)

curve(dnorm(x), add = TRUE, lwd = 2, lty = 2, col = "black") # Standard normal reference

legend("topleft",
   legend = c("Original", "Transformed", "Original Density", "OSKT Density", "Standard Normal"),
   col = c(rgb(0.2,0.4,0.8,0.6), rgb(0.8,0.3,0.3,0.6), "blue", "red", "black"),
   lwd = c(10, 10, 2, 2, 2), lty = c(1, 1, 1, 1, 2), bty = "n")

X_mean <- mean(x_orig)
X_sd   <- sd(x_orig)

res_back <- backosktfast(
              Z = x_transformed,
              X_mean = X_mean, X_sd = X_sd,
              g = g_star, h = h_star,
              method = "brent")

x_recovered <- res_back$X_orig
head(x_recovered, 5)

breaks <- pretty(range(c(x_orig, x_transformed, x_recovered)), n = 30)
hist(x_orig, breaks = breaks, freq = FALSE, col = rgb(0.2, 0.4, 0.9, 0.4),
  border = "white", main="OSKT Transformation & Back Transformation", xlab="Value")
hist(x_transformed, breaks = breaks, freq = FALSE, col = rgb(0.8, 0.3, 0.3, 0.4), 
  border = "white", add=TRUE)
hist(x_recovered, breaks = breaks, freq = FALSE, col = rgb(0.2,0.8,0.2,0.4), 
  border = "white", add=TRUE)

legend("topleft", legend = c("Original","Transformed","Back-transformed"),
       fill = c(rgb(0.2,0.4,0.8,0.4), rgb(0.8,0.3,0.3,0.4), rgb(0.2,0.8,0.2,0.4)))
       
(all.equal(x_orig, x_recovered, tolerance = 1e-6)) # Comparison to the originals with a tolerance


ok <- is.finite(x_orig) & is.finite(x_recovered) # Remove any non-finite values

xo <- x_orig[ok]
xr <- x_recovered[ok]
err <- xr - xo #Error

# Performance metrics
MAE  <- mean(abs(err))
MAXE <- max(abs(err))
MSE  <- mean(err^2)
RMSE <- sqrt(MSE)

# Correlation and R^2
COR  <- cor(xo, xr)
R2   <- COR^2

# Linear fit (recovered ~ original)
fit <- lm(xr ~ xo)
slope     <- coef(fit)[2]
intercept <- coef(fit)[1]

# Summary table
back_stats <- data.frame(
  MSE = MSE, RMSE = RMSE, MAE = MAE, MaxError = MAXE,
  Correlation= COR, R2 = R2, Slope = slope, Intercept  = intercept)

round(t(back_stats), 8)

# Regression line for diagnostic
plot(xo, xr,
     pch = 16, cex = 0.6, col = rgb(0.2, 0.2, 0.7, 0.4),
     xlab = "Original values", ylab = "Back-transformed values",
     main = "OSKT Back-transformation Accuracy")

abline(0, 1, col = "red", lwd = 2, lty = 3) # 1:1 reference line
abline(fit, col = "darkgreen", lwd = 3) # Fitted regression line

legend("topleft",
       legend = c("y = x (ideal)", "Fitted line"),
       col = c("red", "darkgreen"),
       lwd = 2, lty = c(2,1), bty = "n")

# Generate left-skewed data
set.seed(12)
x_orig <- groupcompare::ghdist(n=300, A=0, B=1, g=-0.49, h=0)
x_transformed <- osktfast(x_orig)$transformed

breaks <- pretty(range(c(x_orig, x_transformed)), n = 25)
h_orig  <- hist(x_orig, breaks = breaks, plot = FALSE)
h_trans <- hist(x_transformed, breaks = breaks, plot = FALSE)

d_orig  <- density(x_orig); d_trans <- density(x_transformed)

ymax <- max(c(h_orig$density, h_trans$density, d_orig$y,d_trans$y, dnorm(0)))
hist(x_orig, breaks = breaks, freq = FALSE, ylim = c(0, ymax * 1.05), 
     col = rgb(0.2, 0.4, 0.8, 0.4), border = "white", 
     main = "Before and After OSKT Transformation", xlab = "Value")
lines(d_orig, col = "blue", lwd = 2)

hist(x_transformed, breaks = breaks, freq = FALSE,
     col = rgb(0.8, 0.3, 0.3, 0.4), border = "white", add = TRUE)
lines(d_trans, col = "red", lwd = 2)

curve(dnorm(x), add = TRUE, lwd = 2, lty = 2, col = "black") # Standard normal reference

legend("topleft",
   legend = c("Original", "Transformed", "Original Density", "OSKT Density", "Standard Normal"),
   col = c(rgb(0.2,0.4,0.8,0.6), rgb(0.8,0.3,0.3,0.6), "blue", "red", "black"),
   lwd = c(10, 10, 2, 2, 2), lty = c(1, 1, 1, 1, 2), bty = "n")

# Normality tests with BC, YJ and OSKT 
x_bc <- boxcox(x_orig, makepositive=TRUE)$transformed # Box-Cox transformation
x_yj <- yeojohnson(x_orig)$transformed  # Yeo-Johnson transformation
x_oskt <- osktfast(x_orig)$transformed # OSKT

# Normality tests and moments
get_stats <- function(x) {
  x <- x[is.finite(x)]
  c(
    Skew = mean((x - mean(x))^3) / sd(x)^3,
    Kurt = mean((x - mean(x))^4) / sd(x)^4 - 3,
    SW  = shapiro.test(x)$p.value,
    ZA  = zatest(x, nsim=100)$p.value,
    CVM = cvmtest(x)$p.value,
    PPM = unname(pearsonp(x)$statistic)
  )
}

# Normality results table
pval_table <- rbind(
  ORG  = get_stats(x_orig),
  BC  = get_stats(x_bc),
  YJ   = get_stats(x_yj),
  OSKT = get_stats(x_oskt)
)

pval_table <- as.data.frame(round(pval_table, 4))
pval_table
