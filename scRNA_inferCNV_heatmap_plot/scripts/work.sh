#######################commands below are to reproduce the results ###################################
echo "infercnv final profile"
ln -s ../scRNA_inferCNV_heatmap/output/FEL011046.infercnv.observations.RDS ./
ln -s ../scRNA_inferCNV_heatmap/output/FEL001010.infercnv.observations.RDS ./

echo "merge cell type annotaiton to four broad categories: immune, stromal, cancer and normal epithelial cells"
# output file is FEL011046.cell_metadata.txt
cp ../scRNA_meta/output/FEL001046_scRNA.metadata.clinical.txt FEL001046_scRNA.metadata.clinical.txt
python FEL011046.cell_metadata.py --input FEL001046_scRNA.metadata.clinical.txt
## create a copy for FEL001010
ln -s FEL011046.cell_metadata.txt FEL001010.cell_metadata.txt


echo "plot infercnv heatmap for 35 10x patients and 10 icell8 patients"
# input 1: hg19.RefSeq.NM_pos_unique_sort.txt, gene coordinate from Lance
# OR4F5   1       69091   70008
# OR4F16  1       367659  368597
# input 2: FEL0*.cell_metadata.txt, cell annotation
# input 3: FEL0*.infercnv.observations.RDS, infercnv profile (final profile, not hmm or preliminary)
# output is FEL0*_CNA_infer_NoClust_Celltype.png, which is raw figure for Figure 1c and Supplementary Figure 1b.
sbatch infercnv_heatmap_individual.sh

