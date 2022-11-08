#!/bin/tcsh
#BSUB -W 100
#BSUB -n 1
#BSUB -o XBplot.%J
#BSUB -e XBplot.%J

conda activate /usr/local/usrapps/jcdietri/conda/env_XB
python plot_results_wMap_IRENE.py
conda deactivate

cp -r plots/ ../WEBSITE_public/IRENE2_plots/
cd ../WEBSITE_public/

git init
git branch -M main
git add -A
git commit -am "Irene2_plots"
git push
