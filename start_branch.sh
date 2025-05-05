#!bash

cur=$(git symbolic-ref --short HEAD)

git checkout $1 2>/dev/null || git checkout -b $1


Rscript -e 'library(targets);tar_config_set(script = "_targets.R", store = "_targets_'$1'", project = "'$1'")'

sed -i '' 's/TAR_PROJECT = "'$cur'"/TAR_PROJECT = "'$1'"/g' report.qmd

if [[ "$cur" == "main" ]]; then
  par_tar='_targets'
else
  par_tar='_targets_'$cur
fi


rsync -av -f"+ */" -f"- *" $par_tar/ _targets_$1/
ln -s $par_tar/objects/* _targets_$1/objects/
ln -s $par_tar/meta/* _targets_$1/meta/
ln -s $par_tar/workspaces/* _targets_$1/workspaces/


