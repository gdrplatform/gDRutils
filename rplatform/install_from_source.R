install.packages("git2r", repos = paste0("https://mran.microsoft.com/snapshot/", Sys.getenv("MRAN_SNAPSHOT_DATE")))

## Uncomment following code to install package(s) directly from Bitbucket repository
## SSH keys should be copied to container before installation
ssh_keys <- git2r::cred_ssh_key(file.path("/home/rstudio/.ssh/id_rsa.pub"), file.path("/home/rstudio/.ssh/id_rsa"))
.wd <- "/mnt/vol"
.deps <- yaml::read_yaml(file.path(.wd, "rplatform", "DESCRIPTION_dependencies.yaml"))
pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))$pkgs


for (nm in names(pkgs)) {
  rp::installAndVerify(
    install = devtools::install_git,
    url = pkgs[[nm]]$url,
    ref = pkgs[[nm]]$ref,
    credentials = ssh_keys,
    package = nm,
    # version requirement is taken from DESCRIPTION if not specified manually in yaml
    requirement = if (!is.null(pkgs[[nm]][["ver"]])) pkgs[[nm]][["ver"]] else .deps[[nm]],
    subdir = pkgs[[nm]]$subdir,
    upgrade = "never"
  )
}
