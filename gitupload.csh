#!/bin/tcsh
#BSUB -W 100
#BSUB -n 1
#BSUB -o GITsub.%J
#BSUB -e GITsub.%J

git init
git pull

set name = IRENE_adv99
cp -r ../FIRST_DRYRUN5/plots/ $name

git init
git branch -M main
git add -A
git commit -am "$name"
git push
