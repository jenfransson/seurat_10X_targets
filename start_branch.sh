#!bash

cur=$(git symbolic-ref --short HEAD)

git checkout $1 2>/dev/null || git checkout -b $1


Rscript -e 'library(targets);tar_config_set(script = "_targets.R", store = "_targets_'$1'", project = "'$1'")'

sed -i '' 's/TAR_PROJECT = "'$cur'"/TAR_PROJECT = "'$1'"/g' report.qmd



