
## set cores
if (!require(parallel)) install.packages("parallel")

options(Ncpus = parallel::detectCores())

## install pkgs required for 'rp' package
install.packages(c("git2r", "desc", "devtools", "stringr", "withr", "DelayedMatrixStats"))

devtools::install_git(
  url = "https://bitbucket.roche.com/stash/scm/rp/rp-package.git",
  ref = "master",
  upgrade = "never"
)
