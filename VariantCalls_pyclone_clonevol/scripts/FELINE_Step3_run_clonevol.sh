export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/lib/R/library

for ccf in `cat FELINE_clonevol_ccf.list | grep -v "^#"`; do
    echo $ccf
    prefix=${ccf%.txt}
    a=($(echo $ccf | tr '_' "\n"))
    patient=$a
    analysis='ccf'
    echo $prefix
    echo $patient
    echo $analysis
    #if [ $patient == 'FEL044' ];then
    if true ; then

        if [ $patient == 'FEL044' ] || [ $patient == 'FEL045' ] || [ $patient == 'FEL046' ] || [ $patient == 'FEL036' ]; then
            echo "S vs M"
            /home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/Rscript FELINE_clonevol.R $prefix $patient\_S,$patient\_M $analysis 
        else
            echo "S vs E"
            /home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/Rscript FELINE_clonevol.R $prefix $patient\_S,$patient\_E $analysis
        fi
    fi
done
