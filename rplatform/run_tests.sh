#/bin/sh

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>>>>> RUNNING LINTER"
Rscript -e "gDRstyle::lintPkgDirs('/mnt/vol/gDRutils')"

echo ">>>>> RUNNING UNIT TESTS"
Rscript -e "testthat::test_local(path = '/mnt/vol/gDRutils', stop_on_failure = TRUE)"

echo ">>>>> RUNNING DEPENDENCIES CHECK"
Rscript -e "gDRstyle::checkDependencies('/mnt/vol/gDRutils')"

# TODO: fix the issue with R CMD CHECK
#echo ">>>>> RUNNING DEVTOOLS::CHECK()"
#sudo R CMD check --no-build-vignettes --no-manual --no-tests /mnt/vol/gDRwrapper

echo ">>>>>>>> RUNNING CHECK DEPENDENCIES"
Rscript -e "gDRstyle::checkDependencies(dep_path='/mnt/vol/dependencies.yaml', desc_path='/mnt/vol/gDRutils/DESCRIPTION')"
