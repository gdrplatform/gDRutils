# RP package template version >= 0.0.78

## Uncomment following code to install package from source directory
# rp::installAndVerify(
#   install = devtools::install,
#   requirement = sprintf("== %s", desc::desc_get_version("/mnt/vol/package_source_dir")),
#   package = "/mnt/vol/package_source_dir"
# )

install.packages("git2r", repos = paste0("https://mran.microsoft.com/snapshot/", Sys.getenv("MRAN_SNAPSHOT_DATE")))

## Uncomment following code to install package(s) directly from Bitbucket repository
## SSH keys should be copied to container before installation
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))
.wd <- "/mnt/vol"
