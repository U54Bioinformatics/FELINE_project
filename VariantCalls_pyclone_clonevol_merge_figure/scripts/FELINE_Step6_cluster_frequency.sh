export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/lib/R/library
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library/:$R_LIBS

for ccf in `cat FELINE_clonevol_ccf.list | grep -v "^#"`; do
    echo $ccf
    prefix=${ccf%.txt}
    a=($(echo $ccf | tr '_' "\n"))
    patient=$a
    echo $prefix
    echo $patient
    if [ ! -e $prefix\.cluster_frequency.txt ]; then
        /home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/Rscript FELINE_Step6_cluster_frequency.R $prefix $patient
    fi
done

head -n 1 FEL046_pyclone_analysis_sequenza_20_5_0.05_1_10000_none.clonevol_ccf.cluster_frequency.txt > FELINE_WES_cluster_frequency.txt
cat *.cluster_frequency.txt | grep -v "Patient" >> FELINE_WES_cluster_frequency.txt 
