# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#     Jitter analysis
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Example of reading in different starting values of a program in automated way
# and refitting for each starting value

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
tpl_name <- "surp_prod_jitter"
# setwd("convergence/jitter_test") # set the directory where the .tpl and .dat files are

# compile ADMB
compile_admb(fn = tpl_name, verbose = TRUE)

# set initial values and source from external files
cat("-0.6", file = "log_r.dat", sep = "\n")
cat("-3", file = "log_q.dat", sep = "\n")
cat("8", file = "log_K.dat", sep = "\n")
cat("-1", file = "log_sd_cpue.dat", sep = "\n")

# run ADMB
run_admb(fn = tpl_name, verbose = TRUE)

# get parameter estimates (used for jittering)
if (file.exists("estpars.dat")) {
  dat <- read.table("estpars.dat")
  colnames(dat) <- c("inlog_r", "log_r", "inlog_q", "log_q", "inlog_K", "log_K", "inlog_sd_cpue", "log_sd_cpue", "objn")
}

# Delete any existing version of estpars.dat
if (file.exists("estpars.dat")) file.remove("estpars.dat")

# Create header for file so we know the variables.
# sep ["\n" needed for line feed]
cat("inlog_r log_r inlogq log_q inlogK log_K inlog_sd_cpue log_sd_cpue objn", file = "estpars.dat", sep = "\n")

# Define a set of starting values
nrun <- 50 # number of reruns with new values
st_log_r <- dat$log_r + rnorm(nrun, sd = 0.1)
st_log_q <- dat$log_q + rnorm(nrun, sd = 0.1)
st_log_K <- dat$log_K + rnorm(nrun, sd = 0.1)
st_log_sd_cpue <- dat$log_sd_cpue + rnorm(nrun, sd = 0.1)

# Write out each value of the parameters and run ADMB program for each in loop
for (i in 1:length(st_log_r)) {
  cat(st_log_r[i], file = "log_r.dat", sep = "") # write one st value to file
  cat(st_log_q[i], file = "log_q.dat", sep = "") # write one st value to file
  cat(st_log_K[i], file = "log_K.dat", sep = "") # write one st value to file
  cat(st_log_sd_cpue[i], file = "log_sd_cpue.dat", sep = "") # write one st value to file
  if(Sys.info()["sysname"] == "Windows") { # windows
    system(paste0(tpl_name, ".exe")) 
  } else { # most Mac's and Linux
    system(paste0("./", tpl_name)) 
  }
}

# read in and print results to console
jit_res <- read.table("estpars.dat", header = T)
jit_res

# boxplots - are there any weird shapes/outliers?
boxplot(jit_res[, c(2, 4, 6, 8)])

# clean extra files
clean_admb(fn = tpl_name)
if (file.exists("estpars.dat")) file.remove("estpars.dat")
if (file.exists("out.rdat")) file.remove("out.rdat")
if (file.exists("hessian.rdat")) file.remove("hessian.rdat")
if (file.exists("log_K.dat")) file.remove("log_K.dat")
if (file.exists("log_r.dat")) file.remove("log_r.dat")
if (file.exists("log_q.dat")) file.remove("log_q.dat")
if (file.exists("log_sd_cpue.dat")) file.remove("log_sd_cpue.dat")
