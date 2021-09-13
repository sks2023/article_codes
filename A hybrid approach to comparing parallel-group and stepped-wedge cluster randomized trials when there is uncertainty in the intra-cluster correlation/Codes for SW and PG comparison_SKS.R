##### figures_and_tables.R #####################################################
# Author: Samuel Sarkodie (S.K.Sarkodie2@newcastle.ac.uk)                      #
# Last modified: 12 Sept 2021                                                  #
################################################################################

################################################################################
# INSTALL AND LOAD REQUIRED PACKAGES                                           #
################################################################################

#install.packages("pracma")
#install.packages("tidyverse)
#install.packages("truncnorm")
#install.packages("patchwork")
#install.packages("ggplot2)
library(pracma)
library(tidyverse)
library(truncnorm)
library(ggplot2)
library(patchwork)

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR AN INDIVIDUAL RANDOMISED #
# CONTROLLED TRIAL (PRIOR ONLY ON EFFECT)                                      #
################################################################################

# Function to compute n
n_z           <- function(alpha, beta, sigma, delta) {
  2*(qnorm(1 - alpha) + qnorm(1 - beta))^2*(sigma^2)/(delta^2)
}

# Phi Component
phi_component <- function(delta, n, beta, sigma, alpha) {
  pnorm(delta*sqrt(n/(2*sigma^2)) - qnorm(1 - alpha))
}

# psi component
psi_component <- function(delta, theta_mean, theta_sd) {
  dnorm(delta, theta_mean, theta_sd)
}

# Forming the integrand
integrand     <- function(delta, sigma, beta, alpha, n, theta_mean, theta_sd) {
  phi_component(delta, n, beta, sigma, alpha)*
    psi_component(delta, theta_mean, theta_sd)
}

# Integrating the integrand
pos           <- function(n, sigma, beta, alpha, theta_mean, theta_sd,
                          delta_MCID) {
  integrate(integrand, lower = delta_MCID, upper = Inf,
            sigma = sigma, beta = beta, alpha = alpha, n = n,
            theta_mean = theta_mean, theta_sd = theta_sd)
}

# Testing the components of the function
n_z(alpha = 0.025, beta = 0.1, sigma = 1, delta = 0.35) # ~172
phi_component(delta = 0.35, n = 172, beta = 0.1, sigma = 1,
              alpha = 0.025) # ~0.9
psi_component(delta = 0.35, theta_mean = 0.35, theta_sd = 0.01) # Big
psi_component(delta = 0.35, theta_mean = 0.3, theta_sd = 0.01) # Small
integrand(delta = 0.35, sigma = 1, beta = 0.1, alpha = 0.025, n = 172,
          theta_mean = 0.35, theta_sd = 0.01)
integrand(delta = seq(0, 1, 0.1), sigma = 1, beta = 0.1, alpha = 0.025, n = 172,
          theta_mean = 0.35, theta_sd = 1)
pos(n = 172, sigma = 1, beta = 0.1, alpha = 0.025, theta_mean = 0.35,
    theta_sd = 0.01, delta_MCID = 0.35) # ~0.46

# Finding the minimal n to achieve a desired value of PoS
# Using Loops as a first method
nmax          <- 10000 # Useful in case you can't reach the desired level!
desired_pos   <- 0.48
for (n in 1:nmax) {
  if (pos(n, sigma = 1, beta = 0.1, alpha = 0.025, theta_mean = 0.35,
          theta_sd = 0.01, delta_MCID = 0.35)$value >= desired_pos) {
    minimal_n <- n
    break
  }
  message(n)
}

#215 samples needed to achieve a PoS of 0.48

# The alternative is to use uniroot
wrapper       <- function(n, sigma, beta, alpha, theta_mean, theta_sd,
                          delta_MCID, desired_pos) {
  pos(n, sigma, beta, alpha, theta_mean, theta_sd, delta_MCID)$value -
    desired_pos
}
uniroot(wrapper, c(1e-6, 1e6), sigma = 1, beta = 0.1, alpha = 0.025,
        theta_mean = 0.35, theta_sd = 0.01, delta_MCID = 0.35,
        desired_pos = 0.48)$root # Matches the above

# This treats n as continuous (which is fine in reality)
# Does it handle not being able to achieve the desired POS better or worse?
uniroot(wrapper, c(1e-6, 1e6), sigma = 1, beta = 0.1, alpha = 0.025,
        theta_mean = 0.35, theta_sd = 0.1, delta_MCID = 0.35,
        desired_pos = 0.6)$root # worse

# Using system interval of the machine - marginally slower than using
# c(1e-6, 1e6), but increases chance you find the required n if it exists
# when desired PoS = 0.5
uniroot(wrapper, c(.Machine$double.xmin, .Machine$double.xmax), sigma = 1,
        beta = 0.1, alpha = 0.025, theta_mean = 0.35, theta_sd = 0.01,
        delta_MCID = 0.35, desired_pos = 0.5, maxiter = 10000)$root

################################################################################
# FUNCTIONS TO COMPUTE THE EXPECTED POWER FOR AN INDIVIDUAL RANDOMISED         #
# CONTROLLED TRIAL (PRIOR ONLY ON EFFECT)                                      #
################################################################################

# psi component denominator
denominator           <- function(delta, theta_mean1, theta_sd1) {
  dnorm(delta, theta_mean1, theta_sd1)
}

# Integrate the denominator
integrate_denominator <- function(theta_mean1, theta_sd1, delta_MCID) {
  integrate(denominator, lower = delta_MCID, upper = Inf, theta_sd1 = theta_sd1,
            theta_mean1 = theta_mean1)
}

# psi component numerator
numerator             <- function(delta, theta_mean1, theta_sd1, delta_MCID) {
  ifelse(delta > delta_MCID, 1, 0)*denominator(delta, theta_mean1, theta_sd1)
}

# psi component for EP
psi_component_ep      <- function(delta, theta_mean1, theta_sd1, delta_MCID) {
  numerator(delta, theta_mean1, theta_sd1, delta_MCID)/
    integrate_denominator(theta_mean1, theta_sd1, delta_MCID)$value
}

# Integrand for EP
integrand_ep          <- function(delta, sigma, beta, alpha, n, theta_mean1,
                                  theta_sd1, delta_MCID) {
  phi_component(delta, n, beta, sigma, alpha)*
    psi_component_ep(delta, theta_mean1, theta_sd1, delta_MCID)
}

#Expected power
ep                    <- function(n, sigma, beta, alpha, theta_mean1, theta_sd1,
                                  delta_MCID) {
  integrate(integrand_ep, lower = delta_MCID, upper = Inf, sigma = sigma,
            beta = beta, alpha = alpha, n = n, theta_mean1 = theta_mean1,
            theta_sd1 = theta_sd1, delta_MCID = delta_MCID)
}

# Trying the new components with the same delta_MCID of 0.2 for the POS
integrate_denominator(theta_sd1 = 1, theta_mean1 = 1,
                      delta_MCID = 0.2) # This is bigger
integrate_denominator(theta_sd1 = 1, theta_mean1 = 0,
                      delta_MCID = 0.2) # than this
numerator(delta = 1, theta_mean1 = 1, theta_sd1 = 1, delta_MCID = 0.2) # Non-zero
numerator(delta = 0.5, theta_mean1 = 1, theta_sd1 = 1,
          delta_MCID = 0.2) # Non-zero, but smaller
numerator(delta = 0, theta_mean1 = 1, theta_sd1 = 1, delta_MCID = 0.2) # Zero
integrand_ep(delta = 0.35, sigma = 1, beta = 0.1, alpha = 0.025, n = 172,
             theta_mean1 = 0.35, theta_sd1 = 0.01, delta_MCID = 0.35) # Zero
integrand_ep(delta = 0.35, sigma = 1, beta = 0.1, alpha = 0.025, n = 172,
             theta_mean = 0.35, theta_sd = 1, delta_MCID = 0.2) # Non-zero
integrand_ep(delta = 0.35, sigma = 1, beta = 0.1, alpha = 0.025, n = 172,
             theta_mean = 0.35, theta_sd = 0.1,
             delta_MCID = 0.2) # Non-zero and larger
ep(n = 172, sigma = 1, beta = 0.1, alpha = 0.025, theta_mean = 0.35,
   theta_sd = 0.01, delta_MCID = 0.35) # EP is bigger than PoS using same
                                       # parameters

# Minimal n to achieve a desired value of EP
# Loop method
nmax                  <- 1000 # Useful in case you can't reach the desired
                              # level!
desired_ep            <- 0.48
for (n in 1:nmax) {
  if (ep(n, sigma = 1, beta = 0.1, alpha = 0.025, theta_mean = 0.35,
         theta_sd = 0.01, delta_MCID = 0.35)$value >= desired_ep) {
    minimal_n         <- n
    break
  }
  message(n)
}

# 56 samples needed to achieve an EP of 0.48

# The uniroot alternative
wrapper_ep            <- function(n, sigma, beta, alpha, theta_mean, theta_sd,
                                  delta_MCID, desired_ep) {
  ep(n, sigma, beta, alpha, theta_mean, theta_sd, delta_MCID)$value - desired_ep
}

uniroot(wrapper_ep, c(1e-6, 1e6), sigma = 1, beta = 0.1, alpha = 0.025,
        theta_mean = 0.35, theta_sd = 0.01, delta_MCID = 0.35,
        desired_ep = 0.48)$root

uniroot(wrapper_ep, c(.Machine$double.xmin, .Machine$double.xmax), sigma = 1,
        beta = 0.1, alpha = 0.025, theta_mean = 0.35, theta_sd = 0.01,
        delta_MCID = 0.35, desired_ep = 0.9, maxiter = 10000)$root
# Handles a higher desired_ep better

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR AN INDIVIDUAL RANDOMISED #
# CONTROLLED TRIAL (PRIOR ON EFFECT AND STANDARD DEVIATION)                    #
################################################################################

# Generating the truncated normal prior
psi_component_2 <- function(sig, l_bound, u_bound, theta_mean2, theta_sd2) {
  dtruncnorm(sig, l_bound, u_bound, theta_mean2, theta_sd2)
}

# The double integrand with delta and sig(ma) as priors
d_integrand     <- function(delta, sig, l_bound, u_bound, beta, alpha, n,
                            theta_mean, theta_sd, theta_mean2, theta_sd2) {
  phi_component(delta, n, beta, sig, alpha)*
    psi_component(delta, theta_mean, theta_sd)*
    psi_component_2(sig, l_bound, u_bound, theta_mean2, theta_sd2)
}


# Integrating the integrand
d_pos           <- function(n, l_bound, u_bound, beta, alpha, theta_mean,
                            theta_sd, theta_mean2, theta_sd2, delta_MCID) {
  dblquad(d_integrand, 0, Inf, delta_MCID, Inf, dim = 2,
          tol = .Machine$double.eps^0.5, l_bound = l_bound, u_bound = u_bound,
          beta = beta, alpha = alpha, n = n, theta_mean = theta_mean,
          theta_sd = theta_sd, theta_mean2 = theta_mean2,
          theta_sd2 = theta_sd2)
}

# Tests
psi_component_2(sig = 0.12, l_bound = 0, u_bound = 1, theta_mean2 = 0.35,
                theta_sd2 = 0.1)
d_integrand(0.3, 1, 0, Inf, 0.2, 0.05, 100, 0.3, 1, 1, 1)
d_pos(l_bound = 0, u_bound = Inf, n = 172, beta = 0.1, alpha = 0.025,
      theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1, theta_sd2 = 0.1,
      delta_MCID = 0.1) # High
d_pos(l_bound = 0, u_bound = Inf, n = 172, beta = 0.1, alpha = 0.025,
      theta_mean = 0.35, theta_sd = 1, theta_mean2 = 1, theta_sd2 = 0.1,
      delta_MCID = 0.2) # Lower
system.time(d_pos(l_bound = 0, u_bound = Inf, n = 172, beta = 0.1,
                  alpha = 0.025, theta_mean = 0.35, theta_sd = 0.01,
                  theta_mean2 = 1, theta_sd2 = 0.1,
                  delta_MCID = 0.2)) # Fast enough

# Minimal n to achieve a desired value of the 2d PoS
# Loop method
nmax            <- 1000 # Useful in case you can't reach the desired level!
desired_pos     <- 0.48
for (n in 1:nmax) {
  if (d_pos(n, l_bound = 0, u_bound = Inf, beta = 0.1, alpha = 0.025,
         theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1,
         theta_sd2 = 0.1, delta_MCID = 0.35)>= desired_pos) {
    minimal_n   <- n
    break
  }
  message(n)
}

# The Uniroot option
d_wrapper       <- function(n, l_bound, u_bound, beta, alpha, theta_mean,
                            theta_sd, theta_mean2, theta_sd2, delta_MCID,
                            desired_pos) {
  d_pos(n, l_bound, u_bound, beta, alpha, theta_mean, theta_sd,
        theta_mean2, theta_sd2, delta_MCID) - desired_pos
}
uniroot(d_wrapper, c(1e-6, 1e6), l_bound = 0, u_bound = Inf, beta = 0.1,
        alpha = 0.025, theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1,
        theta_sd2 = 0.1,delta_MCID = 0.35, desired_pos = 0.48)$root
uniroot(d_wrapper, c(1e-6, 1e6), l_bound = 0, u_bound = Inf, beta = 0.1,
        alpha = 0.025, theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1,
        theta_sd2 = 0.1, delta_MCID = 0.35, desired_pos = 0.6)$root
uniroot(d_wrapper, c(.Machine$double.xmin, .Machine$double.xmax),
        l_bound = 0, u_bound = Inf, beta = 0.1, alpha = 0.025,
        theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1,
        theta_sd2 = 0.1,delta_MCID = 0.35, desired_pos = 0.5, maxiter = 10000)$root
# Using the machine limits is now slow

# Okay speed
system.time(uniroot(d_wrapper, c(1e-6, 1e6), l_bound = 0, u_bound = Inf,
                    beta = 0.1, alpha = 0.025, theta_mean = 0.35,
                    theta_sd = 0.01, theta_mean2 = 1, theta_sd2 = 0.1,
                    delta_MCID = 0.35, desired_pos = 0.5)$root)
# Long
system.time(uniroot(d_wrapper, c(.Machine$double.xmin, .Machine$double.xmax),
                    l_bound = 0, u_bound = Inf, beta = 0.1, alpha = 0.025,
                    theta_mean = 0.35, theta_sd = 0.01, theta_mean2 = 1,
                    theta_sd2 = 0.1, delta_MCID = 0.35, desired_pos = 0.5,
                    maxiter = 10000)$root)

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A PARALLEL-GROUP CLUSTER #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON EFFECT AND ICC)                        #
################################################################################

# Phi Component
pg_phi_component  <- function(delta, rho, alpha, sigma, N, C) {
  pnorm(delta*sqrt(C*N/(4*(1 + (N - 1)*rho)*sigma^2)) - qnorm(1 - alpha))
}

# Prior for rho (truncated normal distribution)
psi_component_rho <- function(rho, l_bound, u_bound, theta_mean2, theta_sd2) {
  dtruncnorm(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}

# Integrand
pg_integrand      <- function(delta, rho, l_bound, u_bound, sigma, alpha, N, C,
                              theta_mean, theta_sd, theta_mean2, theta_sd2) {
  pg_phi_component(delta, rho, alpha, sigma, N, C)*
    psi_component(delta, theta_mean, theta_sd)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}

#Trying the components
pg_phi_component(delta=3,rho=0.1,alph=0.025,sigma=7.5,N=11,C=50)#~0.91
psi_component(delta=3,theta_mean=3,theta_sd=1)#~0.40
psi_component_rho(rho=0.1,l_bound=0,u_bound=1,theta_mean2=0.1,theta_sd2=0.01)#~39.89

pg_integrand(delta=3, rho=0.1, l_bound=0, u_bound=1, sigma=7.5, alpha=0.025,
             N=11, C=50, theta_mean=3, theta_sd=1, theta_mean2 = 0.1,
             theta_sd2 = 0.01) #~14.52

# PoS for PG-CRT
pg_pos            <- function(C, N, l_bound, u_bound, alpha, theta_mean,
                              theta_sd, sigma, theta_mean2, theta_sd2, delta_MCID) {
  dblquad(pg_integrand, xa = delta_MCID, xb = Inf, ya = 0, yb = 1, dim = 2,
          tol = .Machine$double.eps^0.5, sigma = sigma,
          l_bound = l_bound, u_bound = u_bound, alpha = alpha, N = N, C = C,
          theta_mean = theta_mean, theta_sd = theta_sd,
          theta_mean2 = theta_mean2, theta_sd2 = theta_sd2)
}


# Using the parameters C = 50, N = 11, theta_mean2 = 0.1, theta_sd2 = 0.01, 
#alpha = 0.025, sigma = 7.5, delta = 3, rho = 0.1

pg_pos(N = 11, C = 50, l_bound=0, u_bound=1, alpha = 0.025, theta_mean = 3,
       theta_sd = 1, sigma = 7.5, theta_mean2 = 0.1, theta_sd2 = 0.01,
       delta_MCID = 3)
# ~0.48 for prior on mu and ICC


# Minimal C to achieve a desired value of pg_pos
# Loop Method
Cmax              <- 1000 # Useful in case you can't reach the desired level!
desired_pos       <- 0.48
for (C in 1:Cmax) {
  if (pg_pos(C, l_bound=0, u_bound=1, theta_sd2 = 0.1, theta_mean2 = 0.05, N = 11,
             alpha = 0.025, theta_mean = 3, theta_sd = 1, sigma = 7.5,
             delta_MCID = 3) >= desired_pos) {
    minimal_C     <- C
    break
  }
  message(C)
}

# ICC(theta_mean2 = 0.1, theta_sd = 0.05);Mu(theta_mean = 3, theta_sd = 1) => C = 48


################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A PARALLEL-GROUP CLUSTER #
# RANDOMISED CONTROLLED TRIAL (PRIOR ONLY ON ICC)                              #
################################################################################

#### CHANGE HERE
pg_integrand_icconly      <- function(rho, delta, l_bound, u_bound, sigma,
                                      alpha, N, C, theta_mean2, theta_sd2) {
  pg_phi_component(delta, rho, alpha, sigma, N, C)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}


pg_integrand_icconly(rho=0.1,delta=3,l_bound=0,u_bound=1,theta_mean2=0.1,
                     theta_sd2=0.01,sigma=7.5,alpha=0.025,N=11, C=50)

pg_pos_icconly       <- function(C, N, l_bound, u_bound, theta_mean2, theta_sd2,
                                 alpha, sigma, delta) {
  integrate(pg_integrand_icconly, 0, 1, l_bound = l_bound, u_bound = u_bound,
            theta_mean2=theta_mean2, theta_sd2=theta_sd2,
            alpha = alpha, sigma = sigma, delta = delta, N = N, C = C)$value
}

pg_pos_icconly(C = 50, N = 11, alpha = 0.025,sigma = 7.5, theta_mean2=0.1,
               theta_sd2=0.01, l_bound=0, u_bound= 1, delta = 3)
# ~.9 which is higher than POS with prior on both Mu and ICC

#when delta_MCID assumes negative values, PoS~Power
pg_pos(N = 11, C = 50, l_bound=0, u_bound=1, alpha = 0.025, theta_mean = 3,
       theta_sd = 0.1, sigma = 7.5, theta_mean2 = 0.1, theta_sd2 = 0.01,
       delta_MCID = -3)#~.9 same as pg_icc_only

##### This converges to PoS with the prior on Mu if you send
##### delta_MCID to -Inf and theta_sd1 to 0. E.g.,
pg_pos(N = 11, C = 50, l_bound=0, u_bound=1, alpha = 0.025, theta_mean = 3,
       theta_sd = 0.1, sigma = 7.5, theta_mean2 = 0.1, theta_sd2 = 0,
       delta_MCID = -Inf)

pg_pos(N = 11, C = 50, l_bound=0, u_bound=1, alpha = 0.025, theta_mean = 3,
       theta_sd = 0.01, sigma = 7.5, theta_mean2 = 0.1, theta_sd2 = 0.01,
       delta_MCID = -3)#Zero when theta_sd2 > 0

# Dropping the prior on mu
for (C in 1:Cmax) {
  if (pg_pos_icconly(C, delta = 3, theta_mean2 = 0.1, theta_sd2 = 0.01, N = 11,
                     alpha = 0.025, sigma = 7.5, u_bound = 1, l_bound = 0) >= desired_pos) {
    minimal_C        <- C
    break
  }
  message(C)
}
# theta_mean2 = 0.1, theta_sd2 = 0.01  => C = 16


################################################################################
# FUNCTIONS TO COMPUTE THE EXPECTED POWER FOR A PARALLEL-GROUP CLUSTER         #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON EFFECT AND ICC)                        #
################################################################################

# EP integrand
pg_integrand_ep      <- function(delta, rho, l_bound, u_bound, sigma, alpha, N, C,
                                 theta_mean1, theta_sd1,theta_mean2, theta_sd2,
                                 delta_MCID) {
  pg_phi_component(delta, rho, alpha, sigma, N, C)*
    psi_component_ep(delta, theta_mean1, theta_sd1, delta_MCID)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}


pg_integrand_ep(delta = 3, rho = 0.1, alpha = 0.025, sigma = 7.5, N= 11, C=50,
                theta_mean1 = 3, theta_sd1 = 1, delta_MCID =2, l_bound = 0, u_bound = 1,
                theta_mean2 = 0.1, theta_sd2 = 0.01)

pg_ep                   <- function(C, alpha, sigma, N, theta_mean1, theta_sd1,
                                    delta_MCID, l_bound, u_bound,
                                    theta_mean2, theta_sd2) {
  dblquad(pg_integrand_ep, xa = delta_MCID, xb = Inf, ya = 0, yb = 1, dim = 2,
          tol = .Machine$double.eps^0.5, sigma = sigma, l_bound=l_bound, u_bound=u_bound,
          theta_mean2=theta_mean2, theta_sd2=theta_sd2, C = C, alpha = alpha, N = N,
          theta_mean1 = theta_mean1, theta_sd1 = theta_sd1, delta_MCID = delta_MCID)
}

# Expected power for PG-CRT with prior on the ICC
pg_ep(N = 11, alpha = 0.025, sigma = 7.5, C = 50, theta_mean1 = 3,
      theta_sd1 = 1, delta_MCID = 1.2, l_bound = 0, u_bound = 1,
      theta_mean2 = 0.1, theta_sd2 = 0.01) # ~0.84

desired_ep             <- 0.9
for (C in 1:Cmax) {
  if (pg_ep(C, theta_mean2 = 0.1, theta_sd2 = 0.05, N = 11, alpha = 0.025,
            l_bound=0, u_bound=1, theta_mean = 3,
            theta_sd = 0.04, sigma = 7.5, delta_MCID = 3) >= desired_ep) {
    minimal_C          <- C
    break
  }
  message(C)
}

# ICC(theta_mean2 = 0.1, theta_sd2 = 0.05);Mu(theta_mean = 3, theta_sd = 0.04) => C = 48
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.05);Mu(theta_mean = 3, theta_sd = 0.25) => C = 43
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.05);Mu(theta_mean = 3, theta_sd = 1) => C = 33
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.05);Mu(theta_mean = 3, theta_sd = 4) => C = 19

# Higher variance on the ICC
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.1);Mu(theta_mean = 3, theta_sd = 0.04) => C = 55
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.1);Mu(theta_mean = 3, theta_sd = 0.25) => C = 50
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.1);Mu(theta_mean = 3, theta_sd = 1) => C = 38
# ICC(theta_mean2 = 0.1, theta_sd2 = 0.1);Mu(theta_mean = 3, theta_sd = 4) => C = 21


################################################################################
# FUNCTIONS TO COMPUTE THE EXPECTED POWER FOR A PARALLEL-GROUP CLUSTER         #
# RANDOMISED CONTROLLED TRIAL (PRIOR ONLY ON EFFECT)                           #
################################################################################

#Expected power without prior on the ICC
pg_integrand_ep_muonly <- function(delta, rho, alpha, sigma, N, C, theta_mean,
                                   theta_sd, delta_MCID) {
  pg_phi_component(delta, rho, alpha, sigma, N, C)*
    psi_component_ep(delta, theta_mean, theta_sd, delta_MCID)
}

pg_ep_muonly           <- function(N, rho, alpha, sigma, C, theta_mean,
                                   theta_sd, delta_MCID) {
  integrate(pg_integrand_ep_muonly, lower = delta_MCID, upper = Inf,
            sigma = sigma, rho = rho, C = C, alpha = alpha, N = N,
            theta_mean = theta_mean, theta_sd = theta_sd,
            delta_MCID = delta_MCID)
}

# Expected power for PG-CRT without prior on the ICC
pg_ep_muonly(N = 11, rho = 0.1, alpha = 0.025, sigma = 7.5, C = 50,
             theta_mean = 3, theta_sd = 1, delta_MCID = 3) # ~0.97

################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A STEPPED-WEDGE CLUSTER  #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON EFFECT AND ICC)                        #
################################################################################

# Phi component for SW-CRT
sw_phi_component <- function(delta, rho, alpha, sigma, n, C, X_unique) {
  S              <- nrow(X_unique)
  Ti             <- ncol(X_unique)
  U              <- sum((C/S)*X_unique)
  W              <- sum(((C/S)*colSums(X_unique))^2)
  V              <- sum((C/S)*rowSums(X_unique)^2)
  c_numerator    <- n*((1 + rho*(n*Ti - 1))*(C*U - W) + n*rho*(U^2 - C*V))
  c_denometor    <- C*sigma^2*(1 - rho)*(1 + rho*(n*Ti - 1))
  pnorm(delta*sqrt(c_numerator/c_denometor) - qnorm(1 - alpha))
}

# Integrand
sw_integrand     <- function(delta, rho, alpha, sigma, n, C, X_unique, u_bound,
                             l_bound, theta_mean, theta_sd, theta_mean2, theta_sd2) {
  sw_phi_component(delta, rho, alpha, sigma, n, C, X_unique)*
    psi_component(delta, theta_mean, theta_sd)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}



# PoS for SW-CRT
sw_pos           <- function(n, l_bound, u_bound, sigma, alpha, C,X_unique,
                             theta_mean, theta_sd,theta_mean2, theta_sd2, delta_MCID) {
  dblquad(sw_integrand, xa = delta_MCID, xb = Inf, ya = 0, yb = 1, dim = 2,
          tol = .Machine$double.eps^0.5, sigma = sigma, X_unique = X_unique,
          u_bound = u_bound, l_bound= l_bound, alpha = alpha, n = n, C = C,
          theta_mean = theta_mean, theta_sd = theta_sd,
          theta_mean2=theta_mean2, theta_sd2=theta_sd2)
}

# Tests
X_unique <- rbind(c(0, 1, 1, 1, 1, 1, 1),
                  c(0, 0, 1, 1, 1, 1, 1),
                  c(0, 0, 0, 1, 1, 1, 1),
                  c(0, 0, 0, 0, 1, 1, 1),
                  c(0, 0, 0, 0, 0, 1, 1),
                  c(0, 0, 0, 0, 0, 0, 1))
sw_phi_component(delta = 0.0282, rho = 0.2, alpha = 0.005, sigma = 0.433, n = 132,
                 C = 30, X_unique = X_unique) # ~.80
sw_integrand(delta = 0.0282, rho = 0.2, sigma = 4.33, alpha = 0.005, n = 132,
             C = 132, X_unique = X_unique, theta_mean = 0.05, theta_sd = 0.01,
             u_bound = 1, l_bound = 0, theta_mean2 = 0.2, theta_sd2 = 0.1) # ~27.48
sw_pos(alpha = 0.005, sigma = 0.433, n = 132, C = 30, X_unique = X_unique,
       l_bound = 0, u_bound = 1, theta_mean = 0.0282, theta_sd = 1,
       theta_mean2=0.2, theta_sd2= 1, delta_MCID = 0.0282)

# Minimal n to achieve a desired value of PoS in SW-CRT
# Loop Method
nmax             <- 1000 # Useful in case you can't reach the desired level!
desired_pos      <- 0.48
for (n in 1:nmax) {
  if (sw_pos(n, alpha = 0.005, sigma = 0.433, C = 30, X_unique = X_unique,
             u_bound = 1, l_bound = 0, theta_mean = 0.028, theta_sd = 0.04,
             theta_mean2=0.2, theta_sd2=0.05, delta_MCID = 0.028) >= desired_pos) {
    minimal_n    <- n
    break
  }
  message(n)
}

# ICC(theta_mean2 = 0.2, theta_sd = 0.05);Mu(theta_mean = 0.028, theta_sd = 0.04) => C = 103


################################################################################
# FUNCTIONS TO COMPUTE THE EXPECTED POWER FOR A STEPPED-WEDGE CLUSTER          #
# RANDOMISED CONTROLLED TRIAL (PRIOR ON EFFECT AND ICC)                        #
################################################################################

# EP Integrand

sw_integrand_ep     <- function(delta, rho, alpha, sigma, n, C, X_unique, u_bound,
                             l_bound, theta_mean1, theta_sd1,
                             theta_mean2, theta_sd2, delta_MCID) {
  sw_phi_component(delta, rho, alpha, sigma, n, C, X_unique)*
    psi_component_ep(delta, theta_mean1, theta_sd1, delta_MCID)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}


sw_ep           <- function(n, sigma, alpha, C, X_unique, u_bound, l_bound,
                            theta_mean1, theta_sd1, theta_mean2, theta_sd2, delta_MCID) {
  dblquad(sw_integrand_ep, xa = delta_MCID, xb = Inf, ya = 0, yb = 1, dim = 2,
          tol = .Machine$double.eps^0.5, sigma = sigma, X_unique = X_unique,
          l_bound = l_bound, u_bound = u_bound, alpha = alpha, n = n, C = C,
          theta_mean1 = theta_mean1, theta_sd1 = theta_sd1,
          theta_mean2=theta_mean2, theta_sd2=theta_sd2, delta_MCID = delta_MCID)
}


sw_ep(alpha = 0.005, sigma = 0.433, n = 132, C = 30, X_unique = X_unique,
      l_bound = 0, u_bound = 1, theta_mean = 0.278, theta_sd = 0.1,
      theta_mean2=0.2, theta_sd2=0.1, delta_MCID = 0.0282)


desired_ep      <- 0.8
for (n in 1:nmax) {
  if (sw_ep(n, alpha = 0.005, sigma = 0.433, C = 30, X_unique = X_unique,
            l_bound = 0, u_bound = 1, theta_mean = 0.028, theta_sd = 0.01,
            theta_mean2=0.2, theta_sd2= 0.05, delta_MCID = 0.028) >= desired_ep) {
    minimal_n   <- n
    break
  }
  message(n)
}

# ICC(theta_mean2 = 0.2, theta_sd2 = 0.05);Mu(theta_mean = 0.028, theta_sd = 0.01) => n = 87
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.05);Mu(theta_mean = 0.028, theta_sd = 0.02) => n = 66
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.05);Mu(theta_mean = 0.028, theta_sd = 0.03) => n = 53
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.05);Mu(theta_mean = 0.028, theta_sd = 0.04) => n = 45

# Higher variance on the ICC
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.1);Mu(theta_mean = 0.028, theta_sd = 0.01) => n = 86
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.1);Mu(theta_mean = 0.028, theta_sd = 0.02) => n = 65
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.1);Mu(theta_mean = 0.028, theta_sd = 0.03) => n = 53
# ICC(theta_mean2 = 0.2, theta_sd2 = 0.1);Mu(theta_mean = 0.028, theta_sd = 0.04) => n = 44


################################################################################
# FUNCTIONS TO COMPUTE THE PROBABILITY OF SUCCESS FOR A STEPPED-WEDGE CLUSTER  #
# RANDOMISED CONTROLLED TRIAL (PRIOR ONLY ON ICC)                              #
################################################################################

# Dropping the prior on mu
sw_integrand_icconly <- function(rho, delta, alpha, sigma, n, C, X_unique,
                                 u_bound, l_bound, theta_mean2, theta_sd2) {
  sw_phi_component(delta, rho, alpha, sigma, n, C, X_unique)*
    psi_component_rho(rho, l_bound, u_bound, theta_mean2, theta_sd2)
}

sw_pos_icconly       <- function(C, delta, alpha, sigma, n, X_unique, u_bound,
                                 l_bound, theta_mean2, theta_sd2) {
  integrate(sw_integrand_icconly, 0, 1,
            u_bound = u_bound, l_bound = l_bound, alpha = alpha,
            X_unique = X_unique, sigma = sigma, delta = delta, n = n,
            C = C, theta_mean2=theta_mean2, theta_sd2=theta_sd2)$value
}

sw_pos_icconly(C = 30, delta = 0.028, alpha = 0.005, sigma = 0.433, n = 132,
               X_unique = X_unique, l_bound = 0, u_bound = 1,
               theta_mean2 = 0.2, theta_sd2 = 0.01) #~ 0.80

desired_pos <-0.8
for (n in 1:nmax) {
  if (sw_pos_icconly(n, delta = 0.028, alpha = 0.005, sigma = 0.433, C = 30,
                     X_unique = X_unique, u_bound = 1, l_bound = 0,
                     theta_mean2 = 0.2, theta_sd2 = 0.01) >=
      desired_pos) {
    minimal_n        <- n
    break
  }
  message(n)
}
# u_bound = 1,   l_bound = 0, theta_mean2 = 0.2, theta_sd2 = 0.01  => n = 133

################################################################################
# PLOTS                                                                        #
################################################################################

#Plot for the truncated normal prior TN (0, 1, 0.1, 0.1)
m1=seq(0.0001, 0.9999, length.out = 10000)
den<-dtruncnorm(m1, a=0, b=1, mean = 0.1, sd=0.1)

#Assumed ICC values from the HTA trial reports
assumedicc<-c(0.1, 0.05,	0.05,	0.1,	0.05,	0.5,	0.047,	0.062,	0.05,	0.04,	0.05,	0.002,
              0.01,	0.03,	0.05,	0.026,	0.01,	0.4,	0.06,	0.02,0.05,	0.01,	0.1,
              0.025,	0.05,	0.11,	0.05,	0.08,	0.05,	0.03,	0.006,	0.05,	0.01)
x=c(1:33)

p<-tibble(m1=m1, den=den)
q<-tibble(assumedicc=assumedicc, x=x)

i=ggplot(p, aes(m1, den))+
  geom_line()+
  theme_bw()+ xlab(expression(paste(italic(m)[1], sep = ""))) +
  labs(y="density", title = "TN (0, 1, 0.1, 0.1) prior")

j=ggplot(q, aes(assumedicc, x))+
  geom_point()+
  theme_bw()+ xlab(expression(paste(italic(m)[1], sep = ""))) +
  labs(y="value", title = "Assumed ICC from HTA trials")

j+i
 
# Comparing POS to EP for fixed n with varying C (PG-CRT and SW-CRT)
C           <- seq(4, 80, 2)
pos_pg      <- ep_pg <- pos_sw <- ep_sw <- numeric(length(C))
i           <- 1
for (c in C) {
  pos_pg[i] <- pg_pos(C = c, N = 50, u_bound = 1, l_bound = 0, alpha = 0.025,
                      theta_mean = 3, theta_sd = 1, theta_mean2 = 0.1, theta_sd2=0.1,
                      sigma = 7.5, delta_MCID = 3)
  ep_pg[i]  <- pg_ep(C = c, N = 50, alpha = 0.025, sigma = 7.5, theta_mean = 3,
                     theta_sd = 1, delta_MCID = 3, u_bound = 1, l_bound = 0,
                     theta_mean2 = 0.1, theta_sd2=0.1)
  pos_sw[i] <- sw_pos(C = c, n = 132, alpha = 0.005, sigma = 0.433,
                      X_unique = X_unique, u_bound = 1, l_bound = 0,
                      theta_mean = 0.0282, theta_sd = 0.01,
                      theta_mean2 = 0.2, theta_sd2=0.2,delta_MCID = 0.0282)
  ep_sw[i]  <- sw_ep(C = c, n = 132, alpha = 0.005, sigma = 0.433,
                     X_unique = X_unique, u_bound = 1, l_bound = 0,
                     theta_mean = 0.0282, theta_sd = 0.01,
                     theta_mean2 = 0.2, theta_sd2=0.2, delta_MCID = 0.0282)
  i         <- i + 1
  message(c)
}
data        <- tibble::tibble(Design = rep(c("PG-CRT", "SW-CRT"),
                                           each = 2*length(C)),
                              Type   = rep(c("PoS", "EP", "PoS", "EP"),
                                           each = length(C)),
                              Value  = c(pos_pg, ep_pg, pos_sw, ep_sw),
                              C      = rep(C, 4))
ggplot(data, aes(C, Value, col = Type)) + geom_line() + facet_wrap(~Design) +
  theme_bw() + theme(legend.position = "bottom", legend.title = element_blank())


# Checking when PG-CRT is better than SW-CRT for different choices of 
#theta mean and sd
theta_mean2         <- seq(0.01, 0.4, length.out = 100)
theta_sd2           <- seq(0.01, 0.5, length.out = 100)
delta_st            <- c(0.1, 0.2, 0.3, 0.4, 0.5)
C_N_Ti              <- rbind(c(10, 30, 3),
                             c(10, 60, 6),
                             c(10, 75, 3),
                             c(10, 150, 6),
                             c(10, 150, 3),
                             c(10, 300, 6),
                             c(25, 30, 3),
                             c(25, 60, 6),
                             c(25, 75, 3),
                             c(25, 150, 6),
                             c(25, 150, 3),
                             c(25, 300, 6),
                             c(50, 30, 3),
                             c(50, 60, 6),
                             c(50, 75, 3),
                             c(50, 150, 6),
                             c(50, 150, 3),
                             c(50, 300, 6),
                             c(100, 30, 3),
                             c(100, 60, 6),
                             c(100, 75, 3),
                             c(100, 150, 6),
                             c(100, 150, 3),
                             c(100, 300, 6)) # C, N, Ti
index               <- 1:nrow(C_N_Ti)
scenarios           <- expand.grid(theta_mean2, theta_sd2, delta_st, index)
pos_pg_icconly      <- pos_sw_icconly  <- numeric(nrow(scenarios))
for (i in 1:nrow(scenarios)) {
  C                 <- C_N_Ti[scenarios[i, 4], 1]
  N                 <- C_N_Ti[scenarios[i, 4], 2]
  Ti                <- C_N_Ti[scenarios[i, 4], 3]
  n                 <- N/Ti
  if (Ti == 3) {
    X_unique <- rbind(c(0, 0, 1),
                      c(0, 1, 1))
  } else {
    X_unique <- rbind(c(0, 0, 0, 0, 0, 1),
                      c(0, 0, 0, 0, 1, 1),
                      c(0, 0, 0, 1, 1, 1),
                      c(0, 0, 1, 1, 1, 1),
                      c(0, 1, 1, 1, 1, 1))
  }
  pos_pg_icconly[i] <- try(pg_pos_icconly(C = C, N = N, theta_mean2 = scenarios[i, 1],
                                          theta_sd2 = scenarios[i, 2], l_bound = 0, u_bound=1,
                                          delta = scenarios[i, 3], alpha = 0.025, sigma = 1))
  pos_sw_icconly[i] <- try(sw_pos_icconly(C = C, n = n, delta = scenarios[i, 3],
                                          alpha = 0.025, sigma = 1,
                                          X_unique = X_unique,
                                          theta_mean2 = scenarios[i, 1],
                                          theta_sd2 = scenarios[i, 2], l_bound = 0, u_bound=1))
 if (i%%1000 == 0) {
   message(i)
 }
}



data                <-
  tibble::tibble(Value       = as.numeric(pos_pg_icconly) - as.numeric(pos_sw_icconly),
                 theta_mean2 = scenarios[, 1],
                 theta_sd2   = scenarios[, 2],
                 delta_st    = factor(paste0("effect size = ", scenarios[, 3])),
                 C           = factor(C_N_Ti[scenarios[, 4], 1]),
                 N           = factor(C_N_Ti[scenarios[, 4], 2]),
                 Ti          = factor(C_N_Ti[scenarios[, 4], 3]),
                 N_Ti        = factor(paste0("N = ", N, ", T = ", Ti),
                                      c("N = 30, T = 3",
                                        "N = 60, T = 6",
                                        "N = 75, T = 3",
                                        "N = 150, T = 6",
                                        "N = 150, T = 3",
                                        "N = 300, T = 6")))


#generating plots for C=10, 25, 50, 100
C = unique(data$C)
for(i in C) {
  message("Filtering for cluster size ", i)
  cluster_data = filter(data, C == i)
  print(cluster_data)
  g = ggplot(cluster_data,
             aes(theta_mean2, theta_sd2, fill = Value, z = Value)) + geom_tile() +
    facet_grid(delta_st~N_Ti) +
    theme_bw() + scale_fill_viridis_c() +
    theme(legend.position = "bottom") + xlab(expression(paste(italic(m)[1], sep = ""))) +
    ylab(expression(paste(italic(s)[1], sep = ""))) + labs(fill = "PoS(PG) - PoS(SW)") +
    geom_contour(colour = "white", breaks = 0)+
    labs(title = paste0("Cluster = ", i))
  print(g)
  readline("Hit return for next plot")
}


#Finding the minimum and maximum values of the mean and sd 
#that makes the PG better than SW in terms of PoS

minmax<-filter(data, Value >0)%>%
  select(Value, theta_mean2,theta_sd2)%>%
  summarise(min(theta_mean2), max(theta_mean2), min(theta_sd2), max(theta_sd2))
 
print(minmax)

#Percentage of the 33 HTA trials for which assumed ICC makes the PG is better
a=q%>%
  filter(assumedicc<=0.105)
count(a)/33

