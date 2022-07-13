#!/bin/tcsh
#BSUB -W 100
#BSUB -n 1
#BSUB -o GITsub.%J
#BSUB -e GITsub.%J

set name = IRENE_adv99
cp -r plots/ ../WEBSITE_public/$name
cd ../WEBSITE_public/

git init
git branch -M main
git add -A
git commit -am "$name"
git push
