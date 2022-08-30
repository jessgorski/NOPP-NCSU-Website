#!/bin/tcsh
#BSUB -W 100
#BSUB -n 1
#BSUB -o GITsub.%J
#BSUB -e GITsub.%J

git init
git pull

set name = 8-30-22

git init
git branch -M main
git add -A
git commit -am "$name"
git push
