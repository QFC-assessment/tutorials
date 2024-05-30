# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#     RTMB example - surplus production model
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# this is an example of a Schaefer model
# you can run all the diagnostics starting in line 80

library(RTMB)

## Data and parameters ####
# setwd("")
# data
dat <- read.table("convergence/RTMB/surp_prod2.dat", header = TRUE)
data <- list()
data$years <- dat$Year
data$n_year <- length(data$years)
data$cat <- dat$Catch
data$cpue <- dat$Index
# parameters
par <- list()
par$log_K <- 8
par$log_r <- -0.7
par$log_q <- -3
par$log_F <- rep(0, data$n_year)
par$log_sd_cat <- -3
par$log_sd_cpue <- -0.7


## Model ####
f <- function(par) {
  getAll(data, par)
  jnll <- 0

  # transform parameters
  K <- exp(log_K)
  r <- exp(log_r)
  q <- exp(log_q)
  F <- exp(log_F)
  sd_cpue <- exp(log_sd_cpue)
  sd_cat <- exp(log_sd_cat)

  # initialize vectors
  bio <- numeric(n_year + 1) # biomass
  cat_hat <- numeric(n_year) # predicted catch
  cpue_hat <- numeric(n_year) # predicted cpue/index
  expl_out <- numeric(n_year) # predicted exploitation rate

  # project the model forward
  bio[1] <- K
  for (t in 1:n_year) {
    expl <- 1 / (1 + F[t])
    bio[t + 1] <- bio[t] + r * bio[t] * (1 - bio[t] / K) - expl * bio[t]
    cat_hat[t] <- expl * bio[t]
    expl_out[t] <- expl
    cpue_hat[t] <- q * bio[t]
  }

  # likelihoods
  jnll <- jnll - sum(dnorm(log(cat), log(cat_hat), sd_cat, TRUE))
  jnll <- jnll - sum(dnorm(log(cpue), log(cpue_hat), sd_cpue, TRUE))

  # report out
  REPORT(bio)
  REPORT(cat_hat)
  REPORT(expl_out)
  REPORT(cpue_hat)

  jnll
}

# create objective function
# fix sd_cat = 0.05
map <- list(log_sd_cat = factor(NA))
obj <- MakeADFun(f, par, map = map)
f(par) # run as a normal R function, does it work?
obj$fn() # should provide a value
obj$gr() # are there any NAs/NaNs/0's?

# run model
opt <- nlminb(obj$par, obj$fn, obj$gr, control = list(eval.max = 1e3, iter.max = 1e3))

# check convergence message
opt$convergence
opt$message
# check standard errors
sdr <- sdreport(obj)
sdr


## Jitter analysis ####
doone <- function() {
  fit <- nlminb(obj$par + rnorm(length(obj$par), sd = .1),
    obj$fn, obj$gr,
    control = list(eval.max = 1000, iter.max = 1000)
  )
  c(fit$par, "convergence" = fit$convergence)
}
set.seed(123456)
jit <- replicate(50, doone()) 
# warning message may pop out
boxplot(t(jit[c(1:3, 28, 29), ])) # check leading parameters
boxplot(t(jit[4:27, ])) # check exploitation rates


## Likelihood profile ####
# profile using log_r (change with 'name' argument, in the order of obj$par)
pro <- TMB:::tmbprofile(obj, name = 2)
plot(pro)
confint(pro)


## Fit vs data ####
pl <- obj$report(opt$par)
# catch
plot(data$cat, pch = 19, col = "steelblue")
lines(pl$cat_hat)
# index
plot(data$cpue, pch = 19, col = "steelblue")
lines(pl$cpue_hat)


## Things to think about ####
# way to remove the warning message in jitter test?
# estimating sd of catch?

