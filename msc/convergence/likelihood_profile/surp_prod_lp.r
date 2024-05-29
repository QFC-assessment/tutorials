# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#     Likelihood profile
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# this is an example of likelihood profiling in ADMB
# note that in the tpl, the parameter you want to profile over is fixed (in this case, log_r)
# you'll have to manually change which parameter is fixed to profile across

### Need to add admb to the Environment:
# Linux Mint (probably most other Linux):
#   add this line to end of .profile in Home directory:
#   export PATH=~/admb:$PATH
#   
# Mac:
#   add this line right before "export PATH" at end of .zprofile in Home directory:
#   PATH="~/admb:${PATH}"
# 
# Windows:
#   In System Environment variables:
#   Add C:\ADMB-13.2\bin to Path

# source R code (from R2admb) to run ADMB
source("add_functions/base_funs.r")
source("add_functions/clean_admb.r")
tpl_name <- "surp_prod_lp"
# setwd("convergence/likelihood_profile") # set the directory where the .tpl and .dat files are

# compile ADMB
compile_admb(fn = tpl_name, verbose = TRUE)

# set initial values and source from external files
cat("-0.994883", file = "log_r.dat", sep = "\n")
# cat("-7.73871", file = "log_q.dat", sep = "\n")
# cat("7.94027", file = "log_K.dat", sep = "\n")
# cat("-2.09553", file = "log_sd_cpue.dat", sep = "\n")

# run ADMB
run_admb(fn = tpl_name, verbose = TRUE)

if (file.exists("estpars.dat")) {
  dat <- read.table("estpars.dat")
  colnames(dat) <- c("inlog_r", "log_r", "log_q", "log_K", "log_sd_cpue", "objn")
}

# First delete any existing version of estpars.dat
if (file.exists("estpars.dat")) file.remove("estpars.dat")

# create header for file so we know the variables.
# sep ["\n" needed for line feed]
cat("inlog_r log_r log_q log_K objn", file = "estpars.dat", sep = "\n")

# Define a set of starting values - range around "mean" parameter value
st_log_r <- seq(dat$log_r - 0.5, dat$log_r + 0.5, 0.05)

# Write out each value of log_r and run ADMB program for each in loop
for (i in 1:length(st_log_r)) {
  cat(st_log_r[i], file = "log_r.dat", sep = "") # write one st value to file
  if(Sys.info()["sysname"] == "Windows") { # windows
    system(paste0(tpl_name, ".exe")) 
  } else { # most Mac's and Linux
    system(paste0("./", tpl_name)) 
  }
}

# read in and print results to console
lp_res <- read.table("estpars.dat", header = T)
lp_res

# look at profile across fixed parameter values
plot(lp_res$log_r, lp_res$objn, 
     type = "l", ylab = "NLL", xlab = "log_r")
abline(v = -0.994883, lty = 2) # value estimated from base run

# clean extra files
clean_admb(fn = tpl_name)
if (file.exists("estpars.dat")) file.remove("estpars.dat")
if (file.exists("out.rdat")) file.remove("out.rdat")
if (file.exists("hessian.rdat")) file.remove("hessian.rdat")
if (file.exists("log_K.dat")) file.remove("log_K.dat")
if (file.exists("log_r.dat")) file.remove("log_r.dat")
if (file.exists("log_q.dat")) file.remove("log_q.dat")
if (file.exists("log_sd_cpue.dat")) file.remove("log_sd_cpue.dat")
