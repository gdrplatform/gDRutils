
## set cores
if (!require(parallel)) install.packages("parallel")

options(Ncpus = parallel::detectCores())

## install pkgs required for 'rp' package 
### new version of git2r causes errors with installing by ssh
MRAN_SNAPSHOT_DATE <- Sys.getenv("MRAN_SNAPSHOT_DATE")
install.packages("git2r", repos = paste0("https://mran.microsoft.com/snapshot/", MRAN_SNAPSHOT_DATE))
install.packages(c("desc", "devtools", "stringr", "withr", "DelayedMatrixStats"))

## install rp directly from bitbucket
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))

devtools::install_git(
  url = "ssh://git@bitbucket.roche.com:7999/rp/rp-package.git",
  ref = "master",
  credentials = ssh_keys
)
