clean_admb <- function(fn, which = c("sys", "output")) {
  if (length(which) == 1) {
    if (which == "none") {
      return()
    }
    if (which == "all") which <- c("sys", "input", "output", "gen")
  }

  sys.ext <- c(
    "bar", "bgs", "cpp", "ecm", "eva", "htp", "luu", "mc2", "mcm", "o", "rep", "rhes",
    "luu", "mc2", "mcm", "tpl.bak", "out", "cout", "shess"
  )
  sys.files <- paste(fn, sys.ext, sep = ".")
  gen.files <- list.files(pattern = "_gen(\\.tpl)*")
  sys.other <- c(
    "eigv.rpt", "fmin.log", "variance", "sims",
    "hesscheck", "hessian.bin", "dgs2", "diags",
    paste("admodel", c("dep", "hes", "cov"), sep = "."),
    paste("admb2r", c("log"), sep = "."),
    paste("cmpdiff.tmp"),
    paste("gradfil1.tmp"),
    paste("gradfil2.tmp"),
    list.files(pattern = "xx.*.tmp"),
    list.files(pattern = ".*f1b2list.*"),
    list.files(pattern = paste(fn, "\\.[bpr][0-9]+", sep = ""))
  )
  ## FIXME: clean up abandoned 'buffer' files too
  ## f1b2list etc.
  input.ext <- c("pin", "dat")
  input.files <- paste(fn, input.ext, sep = ".")
  output.ext <- c("log", "cor", "std", "par", "psv", "hst", "prf", "mcinfo", "exe", "obj")
  output.files <- paste(fn, output.ext, sep = ".")
  output.files <- c(output.files, list.files(pattern = "\\.plt$"))
  if ("sys" %in% which) unlink(c(sys.files, sys.other))
  if ("input" %in% which) unlink(input.files)
  if ("output" %in% which) unlink(output.files)
  if ("gen" %in% which) unlink(gen.files)
}
