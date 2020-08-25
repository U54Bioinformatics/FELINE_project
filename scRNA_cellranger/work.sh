echo "Organize FELINE scRNA 10x data folder"
sbatch FELINE_cellranger_premRNA.sh

echo "Folder structure"
# 1. Each folder is a run of BETSY preprocess (with cell ranger); For example, COH039 is a project ID we can associate with a folder here or records in trello or Megatron; Patients included in each project are listed below.
#COH039: FEL011
#COH044: FEL012,13,14,15,16
#COH049: FEL017,19,20,21,22
#COH047: FEL023,34
#COH055: FEL025,26,27,28
#COH057: FEL029,30,31,32,33,34,35,36,37,38,39,40
#COH064: FEL041,42,43
#COH073: FEL044,45,46
# 2. Each folder has two subfolders (cell_ranger and qc) and one cell vs. gene count matrix (read_counts.txt)
# 3. cell_ranger folder contains a txt file (cell_ranger/summary.txt), which summarize cell and UMI stastistics of this run and a bam file (cell_ranger/FEL019_S/outs/possorted_genome_bam.bam), which can be used to do variant calling and RNA velocity analysis
# 4. qc folder contains many txt and figure files that BETSY generated
# 5. read_counts.txt is the raw cell vs. gene count/UMI matrix that can be used as start for any analysis.
ls FELINE_cellranger_premRNA/COH0*_cell_ranger/*.txt
ls -d FELINE_cellranger_premRNA/COH0*_cell_ranger/cell_ranger
ls -d FELINE_cellranger_premRNA/COH0*_cell_ranger/qc
