#!/bin/tcsh
#BSUB -W 100
#BSUB -n 1
#BSUB -o XBcreate.%J
#BSUB -e XBcreate.%J

conda activate /usr/local/usrapps/jcdietri/conda/env_XB
python dryruns_ADCIRC.py
conda deactivate

### Transects ==> j loop
set j=0
while($j<4158)
##echo ==> print, must use $ before variable
  if ( -d T$j) then
    echo "T$j"
    cp ../MYCODE.csh T$j/
    cd T$j/
    bsub < MYCODE.csh
    cd ..
  else
    ##echo "no T$j"
  endif
  @ j++
end

sleep 3600 

conda activate /usr/local/usrapps/jcdietri/conda/env_XB
python plot_results_wMap_IRENE.py
conda deactivate

cp -r plots/ ../WEBSITE_public/IRENE3_plots/
cd ../WEBSITE_public/

# change commit to storm name and track number
git init
git branch -M main
git add -A
git commit -am "Irene3_plots"
git push
