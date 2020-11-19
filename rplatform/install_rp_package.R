
## set cores
if (!require(parallel)) install.packages("parallel")

options(Ncpus = parallel::detectCores())

## install pkgs required for 'rp' package
install.packages(c("git2r", "desc", "devtools", "stringr", "withr", "DelayedMatrixStats"))

## install rp directly from bitbucket
ssh_keys <- git2r::cred_ssh_key(file.path("/root", ".ssh", "id_rsa.pub"), file.path("/root", ".ssh", "id_rsa"))

devtools::install_git(
  url = "ssh://git@bitbucket.roche.com:7999/rp/rp-package.git",
  ref = "master",
  credentials = ssh_keys,
  upgrade = "never"
)
