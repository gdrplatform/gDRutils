#/bin/sh
repo_path="$1"

echo "Executing $0"
echo "Environment: ${rp_env}"
echo "Working directory: `pwd`"
echo "Working directory contains: `ls | tr '\n' ' '`"

# exit when any command fails
set -e

echo ">>>>>>>> RUNNING LINTER"
Rscript -e "gDRstyle::lintPkgDirs('$repo_path')"

echo ">>>>> RUNNING RCMDCHECK"
Rscript -e "rcmdcheck::rcmdcheck('$repo_path', error_on = 'warning', args = '--no-manual')"

echo ">>>>>>>> RUNNING CHECK DEPENDENCIES"
Rscript -e "gDRstyle::checkDependencies(dep_path='/mnt/vol/rplatform/dependencies.yaml', desc_path='/mnt/vol/gDRutils/DESCRIPTION')"
