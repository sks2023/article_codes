##### figures_and_tables.R #####################################################
# Author: Samuel Sarkodie (S.K.Sarkodie2@newcastle.ac.uk)                      #
# Last modified: 02 Apr 2022                                                   #
################################################################################

################################################################################
# INSTALL AND LOAD REQUIRED PACKAGES                                           #
################################################################################
#install.packages("pracma")
#install.packages("tidyverse")
#install.packages("truncnorm")
#install.packages("patchwork")
#install.packages("ggplot2")
#install.packages("ggthemes")
#install.packages("RColorBrewer")
#install.packages("colorspace")
#install.packages("gridExtra")
library(pracma)
library(tidyverse)
library(truncnorm)
library(patchwork)
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(colorspace)
library(gridExtra)

################################################################################
# PRIORS                                                                       #
################################################################################

# Truncated Normal (ICC)
tnorm_prior <- function(rho, mean, sd) {
  dtruncnorm(rho, 0, 1, mean, sd)
}

# Beta (ICC)
beta_prior  <- function(rho, shape1, shape2) {
  dbeta(rho, shape1, shape2)
}

# Gamma (total variance)
gamma_prior <- function(sigma, shape, rate) {
  dgamma(sigma, shape, rate)
}

# Test the various components
tnorm_prior(rho = 0.1, mean = 0.1, sd = 0.01)   # ~39.89
beta_prior(rho = 0.1, shape1 = 1, shape2 = 9)   # ~3.874
gamma_prior(sigma = 7.5, shape = 7.5, rate = 1) # ~0.144

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A PARALLEL-GROUP CLUSTER #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON TOTAL VARIANCE AND ICC)                #
################################################################################

# Frequentist Power
frequentist_power_pg <- function(delta, rho, alpha, sigma, N, C) {
  pnorm(delta * sqrt(C * N / (4 * (1 + (N - 1) * rho)
                              * sigma^2)) - qnorm(1 - alpha))
}

# Integrand (Truncated Normal for ICC and Gamma for total variance)
pg_integrand_t       <- function(rho, sigma, alpha, delta, N, C, mean, sd,
                                 shape, rate) {
  frequentist_power_pg(delta, rho, alpha, sigma, N, C) *
    tnorm_prior(rho, mean, sd) * gamma_prior(sigma, shape, rate)
}

# Integrand (Beta for ICC and Gamma for total variance)
pg_integrand_b       <- function(rho, sigma, alpha, delta, N, C, shape1, shape2,
                                 shape, rate) {
  frequentist_power_pg(delta, rho, alpha, sigma, N, C) *
    beta_prior(rho, shape1, shape2) * gamma_prior(sigma, shape, rate)
}

# Test the various components
frequentist_power_pg(delta = 3, rho = 0.1, alpha = 0.025, sigma = 7.5, N = 11,
                     C = 50)                                  # ~0.91
pg_integrand_t(rho = 0.1, sigma = 7.5, alpha = 0.025, delta = 3, N = 11, C = 50,
               mean = 0.1, sd = 0.01, shape = 7.5, rate = 1)  # ~5.25
pg_integrand_b(rho = 0.1, sigma = 7.5, alpha = 0.025, delta = 3, N = 11, C = 50,
               shape1 = 1, shape2 = 9, shape = 7.5, rate = 1) # ~ 0.51

# PoS for PG-CRT (Truncated Normal for ICC and Gamma for total variance)
pg_pos_t            <- function(C, N, alpha, delta, mean, sd, shape, rate) {
  dblquad(pg_integrand_t, xa = 0, xb = 1, ya = 0, yb = Inf, dim = 2,
          tol = .Machine$double.eps^0.5, alpha = alpha, delta = delta, N = N,
          C = C, mean = mean, sd = sd, shape = shape, rate = rate)
}

# PoS for PG-CRT (Beta for ICC and Gamma for total variance)
pg_pos_b            <- function(C, N, alpha, delta, shape1, shape2, shape,
                                rate) {
  dblquad(pg_integrand_b, xa = 0, xb = 1, ya = 0, yb = Inf, dim = 2,
          tol = .Machine$double.eps^0.5, alpha = alpha, delta = delta, N = N,
          C = C, shape1 = shape1, shape2 = shape2, shape = shape, rate = rate)
}

# Test the various components, based on the parameters C = 50, N = 11,
# alpha = 0.025, sigma = 7.5, delta = 3, rho = 0.1 from Surr et al.
pg_pos_t(C = 50, N = 11, alpha = 0.025, delta = 3, mean = 0.1, sd = 0.01,
         shape = 7.5, rate = 1) # ~0.86
pg_pos_b(C = 50, N = 11, alpha = 0.025, delta = 3, shape1 = 1, shape2 = 9,
         shape = 7.5, rate = 1) # ~0.87

# Example of how to determine minimal C to achieve a desired value of PoS (loop
# method)
Cmax                <- 1000 # Used in case you can't reach the desired level
desired_pos         <- 0.9
# Truncated Normal prior on ICC and Gamma prior on the total variance
for (C in 1:Cmax) {
  if (pg_pos_t(C, N = 11, alpha = 0.025, delta = 3, mean = 0.1, sd = 0.01,
               shape = 75, rate = 10) >= desired_pos) {
    minimal_C       <- C
    break
  }
  message(C)
}

#TN(0,1,0.1,0.01)  ; Gamma(75,10) => 49
#TN(0,1,0.1,0.045) ; Gamma(30, 4) => 52
#TN(0,1,0.1,1)     ; Gamma(7.5,1) => 62


# Beta prior on ICC and Gamma prior on the total variance
for (C in 1:Cmax) {
  if (pg_pos_b(C, N = 11, alpha = 0.025, delta = 3, shape1 = 0.8, shape2 = 7.2,
               shape = 75, rate = 10) >= desired_pos) {
    minimal_C       <- C
    break
  }
  message(C)
}

#Be(0.8, 7.2); Gamma(75,10)  => 51
#Be(0.1,0.9) ; Gamma(30, 4)  => 50
#Be(1, 1)    ; Gamma(7.5,1)  => 186 Uniform prior on ICC 
                               #and Gamma prior on the total variance


################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A PARALLEL-GROUP CLUSTER #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON ICC ONLY)                              #
################################################################################

# Integrand
pg_integrand_ticc <- function(rho, alpha, delta, sigma, N, C, mean, sd) {
  frequentist_power_pg(delta, rho, alpha, sigma, N, C) *
    tnorm_prior(rho, mean, sd)
}

pg_integrand_bicc <- function(rho, alpha, delta, sigma, N, C, shape1, shape2) {
  frequentist_power_pg(delta, rho, alpha, sigma, N, C) *
    beta_prior(rho, shape1, shape2)
}

# Test the various components
pg_integrand_ticc(rho = 0.1, alpha = 0.025, delta = 3, sigma = 7.5, N = 11,
                  C = 50, mean = 0.1, sd = 0.01)       # ~36.41
pg_integrand_bicc(rho = 0.1, alpha = 0.025, delta = 3, sigma = 7.5, N = 11,
                  C = 50, shape1 = 0.1, shape2 = 0.01) # ~0.073

# PoS for PG-CRT (Truncated Normal for ICC)
pg_pos_ticc       <- function(C, N, alpha, sigma, delta, mean, sd) {
  integrate(pg_integrand_ticc, 0, 1, C = C, N = N, alpha = alpha, sigma = sigma,
            delta = delta, mean = mean, sd = sd)
}

# PoS for PG-CRT (Beta for ICC)
pg_pos_bicc       <- function(C, N, alpha, sigma, delta, shape1, shape2) {
  integrate(pg_integrand_bicc, 0, 1, C = C, N = N, alpha = alpha, sigma = sigma,
            delta = delta, shape1 = shape1, shape2 = shape2)
}

# Test the various components
pg_pos_ticc(C = 50, N = 11, alpha = 0.025, sigma = 7.5, delta = 3, mean = 0.1,
            sd = 0.01)  # ~0.91
pg_pos_bicc(C = 50, N = 11, alpha = 0.025, sigma = 7.5, delta = 3, shape1 = 1,
            shape2 = 9) # ~0.90

# Sample size with prior on the ICC only
for (C in 1:Cmax) {
  if (pg_pos_ticc(C, N = 11, alpha = 0.025, delta = 3, sigma = 7.5,
                  mean = 0.1, sd = 0.01)$value >= desired_pos) {
    minimal_C       <- C
    break
  }
  message(C)
}

# TN(0, 1, 0.1, 0.01)  = 47
# TN(0, 1, 0.1, 0.05)  = 49
# TN(0, 1, 0.1, 0.1)   = 150

# TN(0, 1, 0.15, 0.01)  = 59
# TN(0, 1, 0.15, 0.087) = 63
# TN(0, 1, 0.15, 0.14)  = 71
# TN(0, 1, 0.15, 0.32)  = 107

# TN(0, 1, 0.05, 0.01)  = 35
# TN(0, 1, 0.05, 0.087) = 46
# TN(0, 1, 0.05, 0.14)  = 57
# TN(0, 1, 0.05, 0.32)  = 96

for (C in 1:Cmax) {
  if (pg_pos_bicc(C, N = 11, alpha = 0.025, delta = 3, sigma = 7.5,
                  shape1 = 0.8, shape2 = 7.2)$value >= desired_pos) {
    minimal_C       <- C
    break
  }
  message(C)
}

# Be(0.8, 7.2) [mean = 0.1, var = 0.01]  = 49  clusters
# Be(0.1, 0.9) [mean = 0.1, var = 0.045] = 47  clusters
# Be(1, 1)     [mean = 0.5, var = 0.08]  = 158 clusters

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A STEPPED-WEDGE CLUSTER  #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON VARIANCE AND ICC)                      #
################################################################################

# Frequentist power for SW-CRT
frequentist_power_sw <- function(delta, rho, alpha, sigma, n, C, X_unique) {
  S           <- nrow(X_unique)
  Ti          <- ncol(X_unique)
  U           <- sum((C/S)*X_unique)
  W           <- sum(((C/S)*colSums(X_unique))^2)
  V           <- sum((C/S)*rowSums(X_unique)^2)
  c_numerator <- n*((1 + rho*(n*Ti - 1))*(C*U - W) + n*rho*(U^2 - C*V))
  c_denometor <- C*sigma^2*(1 - rho)*(1 + rho*(n*Ti - 1))
  pnorm(delta*sqrt(c_numerator/c_denometor) - qnorm(1 - alpha))
}

# Integrand (Truncated Normal for ICC and Gamma for total variance)
sw_integrand_t       <- function(rho, sigma, alpha, delta, n, C, X_unique, mean,
                                 sd, shape, rate) {
  frequentist_power_sw(delta, rho, alpha, sigma, n, C, X_unique) *
    gamma_prior(sigma, shape, rate) * tnorm_prior(rho, mean, sd)
}

# Integrand (Beta for ICC and Gamma for total variance)
sw_integrand_b       <- function(rho, sigma, alpha, delta, n, C, X_unique,
                                 shape1, shape2, shape, rate) {
  frequentist_power_sw(delta, rho, alpha,  sigma, n, C,  X_unique) *
    gamma_prior(sigma, shape, rate) * beta_prior(rho, shape1, shape2)
}

# PoS for SW-CRT (Truncated Normal for ICC and Gamma for total variance)
sw_pos_t             <- function(n, C, alpha, delta, X_unique, mean, sd, shape,
                                 rate) {
  dblquad(sw_integrand_t, xa = 0, xb = 1, ya = 0, yb = Inf, dim = 2,
          tol = .Machine$double.eps^0.5, alpha = alpha, delta = delta, n = n,
          C = C, X_unique = X_unique, mean = mean, sd = sd, shape = shape,
          rate = rate)
}

# PoS for SW-CRT (Beta for ICC and Gamma for total variance)
sw_pos_b             <- function(n, C, alpha, delta, X_unique, shape1, shape2,
                                 shape, rate) {
  dblquad(sw_integrand_b, xa = 0, xb = 1, ya = 0, yb = Inf, dim = 2,
          tol = .Machine$double.eps^0.5, alpha = alpha, delta = delta, n = n,
          C = C, X_unique = X_unique, shape1 = shape1, shape2 = shape2,
          shape = shape, rate = rate)
}

# Test the various components, using parameters from O'Grady et al.
# delta = 0.0282, rho = 0.2, sigma = 0.433, alpha = 0.005, n = 132, C = 30
X_unique             <- rbind(c(0, 1, 1, 1, 1, 1, 1),
                              c(0, 0, 1, 1, 1, 1, 1),
                              c(0, 0, 0, 1, 1, 1, 1),
                              c(0, 0, 0, 0, 1, 1, 1),
                              c(0, 0, 0, 0, 0, 1, 1),
                              c(0, 0, 0, 0, 0, 0, 1))

frequentist_power_sw(delta = 0.0282, rho = 0.2, alpha = 0.005, sigma = 0.433,
                     n = 132, C = 30, X_unique = X_unique)  # ~0.80

sw_integrand_t(rho = 0.2, sigma = 0.433, alpha = 0.005, delta = 0.0282, n = 132,
               C = 30, X_unique = X_unique, mean = 0.2, sd = 0.1, shape = 7.5,
               rate = 17.32)                                # ~8.16

sw_integrand_b(rho = 0.2, sigma = 0.433, alpha = 0.005, delta = 0.0282, n = 132,
               C = 30, X_unique = X_unique, shape1 = 2, shape2 = 8, shape = 7.5,
               rate = 17.32)                                # ~6.04

sw_pos_t(n = 132, C = 30, alpha = 0.005, delta = 0.0282, X_unique = X_unique,
         mean = 0.2, sd = 0.1, shape = 7.5, rate = 17.32)   # ~0.77

sw_pos_b(n = 132, C = 30, alpha = 0.005, delta = 0.0282, X_unique = X_unique,
         shape1 = 2, shape2 = 8, shape = 7.5, rate = 17.32) # 0.77

# Example of how to determine minimal C to achieve a desired value of PoS (loop
# method)
nmax                 <- 1000 # Used in case you can't reach the desired level
desired_pos          <- 0.8
# Truncated Normal prior on ICC and Gamma prior on the total variance
for (n in 1:nmax) {
  if (sw_pos_t(n, C = 30, alpha = 0.005, delta = 0.0282, X_unique = X_unique,
               mean = 0.2, sd = 0.01, shape = 17.32, rate = 40) >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}


#TN(0,1,0.2,0.01)  ; Gamma(17.32,  40) => 139
#TN(0,1,0.2,0.05)  ; Gamma(4.33,   10) => 157
#TN(0,1,0.2,1)     ; Gamma(0.173, 0.4) => 47


# Beta prior on ICC and Gamma prior on the total variance
for (n in 1:nmax) {
  if (sw_pos_b(n, C = 30, alpha = 0.005, delta = 0.0282, X_unique = X_unique,
               shape1 = 3, shape2 = 12, shape = 17.32, rate = 40) >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}

#Be(3,    12)    ; Gamma(17.32,  40) => 139
#Be(0.4, 1.6)    ; Gamma(4.33,   10) => 146

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A STEPPED-WEDGE CLUSTER  #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON ICC ONLY)                              #
################################################################################

sw_integrand_ticc <- function(rho, alpha, sigma, delta, n, C, X_unique, mean,
                              sd) {
  frequentist_power_sw(delta, rho, alpha, sigma, n, C, X_unique) *
    tnorm_prior(rho, mean, sd)

}

sw_integrand_bicc <- function(sigma, rho, alpha, delta, n, C, X_unique, shape1,
                              shape2) {
  frequentist_power_sw(delta, rho, alpha, sigma, n, C, X_unique) *
    beta_prior(rho, shape1, shape2)
}

sw_pos_ticc       <- function(n, C, alpha, sigma, delta, X_unique, mean, sd) {
  integrate(sw_integrand_ticc, 0, 1, n = n, C = C, alpha = alpha, sigma = sigma,
            delta = delta, X_unique = X_unique, mean = mean, sd = sd)
}

sw_pos_bicc       <- function(n, C, alpha, sigma, delta, X_unique, shape1,
                              shape2) {
  integrate(sw_integrand_bicc, 0, 1, n = n, C = C, alpha = alpha, sigma = sigma,
            delta = delta, X_unique = X_unique, shape1 = shape1,
            shape2 = shape2)
}

# Test the various components
sw_pos_ticc(n = 132, C = 30, alpha = 0.005, sigma = 0.433, delta = 0.0282,
            X_unique = X_unique, mean = 0.2, sd = 0.1)   # ~0.81

sw_pos_bicc(n = 132, C = 30, alpha = 0.005, sigma = 0.433, delta = 0.0282,
            X_unique = X_unique, shape1 = 2, shape2 = 8) # ~0.80

# Sensitivity analysis (Correctly specified prior):ICC only 
for (n in 1:nmax) {
  if (sw_pos_bicc(n, C = 30, alpha = 0.005, delta = 0.0278, X_unique = X_unique,
                  shape1 = 3, shape2 = 12, sigma = 0.426)$value >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}

#Be(3, 12)    = 131
#Be(0.4, 1.6) = 125

for (n in 1:nmax) {
  if (sw_pos_ticc(n, C = 30, alpha = 0.005, delta = 0.0278, X_unique = X_unique,
                  mean = 0.2, sd = 0.01, sigma = 0.426)$value >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}

#TN (0, 1, 0.2, 0.01) = 131
#TN (0, 1, 0.2, 0.05) = 131
#TN (0, 1, 0.2, 1)    = 91

#Sensitivity Analysis with prior on ICC only
for (n in 1:nmax) {
  if (sw_pos_ticc(n, C = 30, alpha = 0.005, delta = 0.0278, X_unique = X_unique,
               mean = 0.25, sd = 0.01, sigma = 0.426)$value >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}

#TN(0, 1, 0.25, 0.01)  = 123
#TN(0, 1, 0.25, 0.175) = 119
#TN(0, 1, 0.25, 0.3)   = 109
#TN(0, 1, 0.25, 0.8)   = 92

#TN(0, 1, 0.15, 0.01)  = 139
#TN(0, 1, 0.15, 0.175) = 129
#TN(0, 1, 0.15, 0.3)   = 116
#TN(0, 1, 0.15, 0.8)   = 94
