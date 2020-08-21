export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/lib/R/library
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library/:$R_LIBS

for ccf in `cat FELINE_clonevol_ccf.list | grep -v "^#"`; do
    echo $ccf
    prefix=${ccf%.txt}
    a=($(echo $ccf | tr '_' "\n"))
    patient=$a
    echo $prefix
    echo $patient
    /home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/Rscript FELINE_Step4_run_fishplot.R $prefix $patient
    convert -density 100 -quality 100 $prefix\.clonevol_figure4_fishplot.pdf $prefix\.clonevol_figure4_fishplot.png
done
