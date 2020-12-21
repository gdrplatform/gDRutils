# RP package template version >= 0.0.78

# The base image makes use of the stable GRAN repository corresponding to image's R version as a source of packages.
# If you wish to install custom package versions from other source you'll need to
# modify the "repos" option or provide the repository url explicitly during package installation.

# To install and verify packages use rp::installAndVerify function.
# Version requirements can be defined with "requirement" argument:
# -- install = install.packages
# '*' - any version
# '==0.1' - package version equal to 0.1
# '>=0.1' - package version greater than or equal to 0.1
# -- install = install_github
# 'r:tag_name' - reference github tag
# 's:sha1' - github SHA1

# Install packages using multiple cores
if (!require(parallel)) install.packages("parallel")
options(Ncpus = parallel::detectCores())

.wd <- "/mnt/vol"

# Extract dependencies from DESCRIPTION file
deps <- yaml::read_yaml(file.path(.wd, "rplatform", "DESCRIPTION_dependencies.yaml"))
deps <- deps[!(names(deps) %in% dont.install)]
# packages needed in
# deps <- rbind(deps, data.frame(type = rep("Suggests", 3),
#                               package = c("attempt", "urltools", "config"),
#                               version = rep("*", 3)))

for (nm in names(deps)) {
  devtools::install_version(
    package = nm,
    version = if (deps[[nm]] == "*") NULL else deps[[nm]])
}
