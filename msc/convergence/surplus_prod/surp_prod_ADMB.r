# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#
#     Surplus production model in ADMB
#
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# Schaefer production model for South African hake
# running ADMB through R without any packages

# source R code (from R2admb) to run ADMB
source("add_functions/base_funs.r")
source("add_functions/clean_admb.r")
tpl_name <- "surp_prod"
# setwd("convergence/surplus_prod") # set the directory where the .tpl and .dat files are

# compile ADMB
compile_admb(fn = tpl_name, verbose = TRUE)

# run ADMB
run_admb(fn = tpl_name, verbose = TRUE)

# look at hessian
hess <- dget("hessian.rdat")$hessian
# is the Hessian positive definite?
ifelse(is.nan(max(hess)), 
        "Hessian is not positive definite", 
        "Hessian is positive definite")


# Plot fit vs data
# before looking at this, do jitter test and likelihood profile
res <- dget("out.rdat")
# index
plot(res$obs_cpue, pch = 19)
lines(res$est_cpue, lwd = 2)
# catch
plot(res$obs_cat, pch = 19)
lines(res$est_cat, lwd = 2)


# clean extra files
clean_admb(fn = tpl_name)
if (file.exists("estpars.dat")) file.remove("estpars.dat")
if (file.exists("out.rdat")) file.remove("out.rdat")
if (file.exists("hessian.rdat")) file.remove("hessian.rdat")
if (file.exists("log_K.dat")) file.remove("log_K.dat")
if (file.exists("log_r.dat")) file.remove("log_r.dat")
if (file.exists("log_q.dat")) file.remove("log_q.dat")
if (file.exists("log_sd_cpue.dat")) file.remove("log_sd_cpue.dat")
