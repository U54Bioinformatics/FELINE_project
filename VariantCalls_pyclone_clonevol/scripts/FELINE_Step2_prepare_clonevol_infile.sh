for indir in `cat FELINE_clonevol_prepare_table.list`; do
    a=($(echo $indir | tr '_' "\n"))
    patient=$a
    #if [ $patient == 'FEL025' ];then
    if true ; then
        echo $indir
        ## from Jeff's analysis table: cluster.frequency.final.txt
        #python FELINE_clonevol_prepare_table.py --pyclone_dir $indir --driver Cancer_genes.ID.list
        ## from pyclone table: tables/pyclone.table_cluster.final.txt
        python FELINE_clonevol_prepare_table_from_pyclone_table.py --pyclone_dir $indir --driver Cancer_genes.ID.list 
    fi
done
