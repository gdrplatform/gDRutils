install.packages("git2r", repos = paste0("https://mran.microsoft.com/snapshot/", Sys.getenv("MRAN_SNAPSHOT_DATE")))

## Uncomment following code to install package(s) directly from Bitbucket repository
.wd <- "/mnt/vol"
.deps <- yaml::read_yaml(file.path(.wd, "rplatform", "DESCRIPTION_dependencies.yaml"))
pkgs <- yaml::read_yaml(file.path(.wd, "rplatform", "git_dependencies.yml"))$pkgs


for (nm in names(pkgs)) {
  rp::installAndVerify(
    install = devtools::install_git,
    url = pkgs[[nm]]$url,
    ref = pkgs[[nm]]$ref,
    package = nm,
    # version requirement is taken from DESCRIPTION if not specified manually in yaml
    requirement = if (!is.null(pkgs[[nm]][["ver"]])) pkgs[[nm]][["ver"]] else .deps[[nm]],
    subdir = pkgs[[nm]]$subdir,
    upgrade = "never"
  )
}
