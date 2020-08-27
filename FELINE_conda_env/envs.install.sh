echo "install BETSY"
git clone ssh://jchen2@129.106.148.73:53849/home/changlab/repos/changlab.git ./
python setup.py build
python setup.py install --prefix /home/jichen/software/Python_lib
export PATH=$PATH:/home/jichen/software/Python_lib64/lib64/python2.7/site-packages/genomicode/bin/

wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p ./install
./install/bin/conda install -c conda-forge -c r -c bioconda numpy matplotlib pygraphviz rpy2 xlwt openpyxl xlrd pil pandas parallel
./install/bin/python setup.py install --prefix /home/jichen/software/Python_lib
./install/bin/conda install -c conda-forge -c r -c bioconda samtools bcftools picard vcftools bedtools tophat
./install/bin/conda install -c conda-forge -c r -c bioconda trimmomatic fastqc bowtie bowtie2 bwa star
./install/bin/conda install -c conda-forge -c r -c bioconda gatk platypus-variant varscan snpeff pindel
./install/bin/conda install -c conda-forge -c r -c bioconda rsem htseq rseqc
./install/bin/conda install -c conda-forge -c pkgs/main -c pkgs/free python-dateutil=2.6.1 psutil=5.6.2 pil=1.1.7
export PATH=$PATH:/home/jichen/software/BETSY/install/bin
./install/bin/python setup.py install --prefix /home/jichen/software/Python_lib

#install GATK
cd /home/jichen/software/BETSY/install/opt/gatk-3.8
tar jxf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2
cd GenomeAnalysisTK-3.8-1-0-gf15c1c3ef
mv GenomeAnalysisTK.jar ../
cd ..
rm -R GenomeAnalysisTK-3.8-1-0-gf15c1c3ef

#Variant calling
./install/bin/conda install -c conda-forge -c r -c bioconda picard=2.18.4 manta varscan
./install/bin/conda install -c conda-forge -c r -c bioconda snpeff strelka snp-pileup sambamba samblaster bam-readcount 
./install/bin/conda install -c shahcompbio -c matnguyen -c pwwang -c aroth85 museq strelka1 cnvnator pyclone
./install/bin/conda install -c conda-forge -c r -c bioconda sequenza-utils htslib 
./install/bin/conda install -c conda-forge -c r -c bioconda django-model-utils
wget https://raw.githubusercontent.com/juanpex/django-model-report/master/model_report/arial10.py
mv arial10.py ~/software/Python_lib/lib/python2.7/site-packages/genomicode/
cp ~/software/Python_lib/lib/python2.7/site-packages/genomicode/arial10.py ~/software/BETSY/install/lib/python2.7/site-packages/
#gdbm needed by ssGSEA
./install/bin/conda install python-gdbm
python -c "import gdbm"
#gseabase
./install/bin/conda install -c bioconda bioconductor-gseabase

./install/bin/conda install -c biobuilds r-sequenza
#R libraries are install here
#~/software/BETSY/install/lib/R/library/

#install cluster 3.0 for pyclone
cd /home/jichen/software/BETSY/install/pkgs
wget http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster-1.58.tar.gz
uz cluster-1.58.tar.gz
cd cluster-1.58
./configure --prefix=/home/jichen/software/BETSY/install/ --without-x
make
make install
#install PIL library for pyclone
./install/bin/conda install -c anaconda pillow
cp /home/jichen/software/BETSY/install/lib64/python2.7/site-packages/genomicode/MS\ PGothic.ttf /home/jichen/software/BETSY/install/lib/python2.7/site-packages/genomicode/


#fix sequenza and facets R problems
./install/bin/conda install 'icu=58.*'
./install/bin/conda install -c conda-forge -c r -c bioconda rpy2
#remove system R and use R in conda that is compitable with rpy2
#test in python
#import rpy2.robjects, no Segmentation fault (core dumped)
#sequenza and facets (confirmed) use chromosome numbers (1-22,X,Y,M, not chr1-chr22).
#need remove chr after run step1. 
#sed 's/chr//g' 101024___101025.pileup.txt > ../test_FACETSPileupFolder2/101024___101025.pileup.txt
#sed 's/chr//g' 101024___101026.pileup.txt > ../test_FACETSPileupFolder2/101024___101026.pileup.txt
#sed 's/chr//g' 101024___101027.pileup.txt > ../test_FACETSPileupFolder2/101024___101027.pileup.txt

#install sequenza in anothre env
./install/bin/conda create -n sequenza -c conda-forge -c r -c bioconda -c biobuilds r-sequenza=3.0.0
#it works but running into errors with data
#we switch to version 2.1.2 with Jeff's edit
#./install/bin/conda install conda-build
#./install/bin/conda index ./Sequenza_BETSY_fi
cd Sequenza_BETSY_fix/
~/software/BETSY/install/envs/sequenza_2_1_2/bin/R CMD INSTALL sequenza

#scRNA
./install/bin/conda install -c bioconda bioconductor-biomart

#try to install all scRNA preprocess and QC in a single env, 20190603
#./install/bin/conda create -n scRNA -c r -c bioconda r-data.table r-matrix
#./install/bin/conda install -n scRNA -c r -c bioconda -c conda-forge -c derickl r-dplyr r-seurat umap-learn
#./install/bin/conda install -n scRNA -c conda-forge r-beeswarm
#./install/bin/conda install -c conda-forge umap-learn=0.3.8
#./install/bin/conda install -c pkgs/r r-base=3.5.1
./install/bin/conda init
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
./install/bin/conda create -n scRNA -c r -c bioconda -c pkgs/r r-base=3.5.1
./install/bin/conda install -n scRNA -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr
#install develop seurat to fix cluster orderin problems
#in R
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
#install seurat-wrappers
devtools::install_github('satijalab/seurat-wrappers')
#installvelocyto.R
#./install/bin/conda install -n scRNA -c anaconda hdf5
./install/bin/conda install -n scRNA -c conda-forge -c bioconda -c statiskit libboost-dev bioconductor-pcamethods r-hdf5r
cd install/lib
ln -s ~/software/BETSY/install/envs/scRNA/lib/libboost_* ./
R
library(devtools)
install_github("velocyto-team/velocyto.R")
#loomR
module load hdf5/1.10.3
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
library(loomR)
#conos
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
R
devtools::install_github("hms-dbmi/conos")
#MAST
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
BiocManager::install("MAST")
library("MAST")
##fix umap bugs: address 0xfffffffffffffff7, cause 'memory not mapped'
/home/jichen/software/BETSY/install/envs/scRNA/bin/R
remove.packages('RcppParallel')
./install/bin/conda install -n scRNA -c conda-forge r-rcppparallel
./install/bin/conda install -n scRNA -c conda-forge r-circlize

#scRNA_velocityR in a single env
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_velocityR/lib/R/library
./install/bin/conda create -n scRNA_velocityR -c r -c bioconda -c pkgs/r r-base=3.5.1
./install/bin/conda install -n scRNA_velocityR -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr
#in R, install seurat and wrappers
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
install.packages("devtools")
library(devtools)
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
install.packages("dplyr")
devtools::install_github('satijalab/seurat-wrappers')
quit() 
#exit R: 
#install velocityR package
./install/bin/conda install -n scRNA_velocityR -c conda-forge -c bioconda -c statiskit libboost-dev bioconductor-pcamethods r-hdf5r
./install/bin/conda install -n scRNA_velocityR -c anaconda boost
#in R
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
install_github("velocyto-team/velocyto.R")
quit()
#exit R
#install loomR 
module load hdf5/1.10.3
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
devtools::install_github(repo = "hhoeflin/hdf5r")
devtools::install_github(repo = "mojaveazure/loomR", ref = "develop")
#install r-tidyverse
./install/bin/conda install -n scRNA_velocityR -c r r-tidyverse
#reinstall edited version of velocity.R in /home/jichen/software/Single_cell_RNAseq/RNAvelocity
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_velocityR/lib/R/library
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
install.packages("callr")
install.packages("devtools")
library(devtools)
install_github("JinfengChen/velocyto.R")
#fix umap bugs: address 0xfffffffffffffff7, cause 'memory not mapped'
/home/jichen/software/BETSY/install/envs/scRNA_velocityR/bin/R
remove.packages('RcppParallel')
./install/bin/conda install -n scRNA_velocityR -c conda-forge r-rcppparallel

#scRNA_QC, remove low quality cells, doublets, maybe normalization too.
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_QC/lib/R/library
./install/bin/conda create -n scRNA_QC -c r -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.5.1 r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr bioconductor-scater
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_QC/lib/python3.6/site-packages/
./install/bin/conda install -n scRNA_QC -c grst -c conda-forge python-annoy scrublet

#
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/lib/R/library
./install/bin/conda create -n scRNA_QC_R3.6 -c r -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.6 r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr bioconductor-scater bioconductor-scran bioconductor-biocsingular
/home/jichen/software/BETSY/install/envs/scRNA_QC_R3.6/bin/R
#install the lastest version scater 1.14 in R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scater")
library(scater)
sessionInfo()
#install scds
BiocManager::install("scds", version = "3.10")
#seurat wrappers
install.packages("devtools")
#install dplyr to fix Rcpp errors
install.packages("dplyr", type = "source")
library(dplyr)
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
#install seurat-wrappers
devtools::install_github('satijalab/seurat-wrappers')
./install/bin/conda install -n scRNA_QC_R3.6 -c r -c bioconda -c pkgs/r -c conda-forge umap-learn=0.3.8
./install/bin/conda install -n scRNA_seurat_wrapper -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr bioconductor-batchelor


#velocyto
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/scRNA_velocyto/lib/python3.6/site-packages/
./install/bin/conda env create -f velocyto.yml -n scRNA_velocyto
#scVelo in the same env
/home/jichen/software/BETSY/install/envs/scRNA_velocyto/bin/pip install -U scvelo
#igragh in the same env
./install/bin/conda install -n scRNA_velocyto -c conda-forge python-igraph
#
./install/bin/conda install -n scRNA_velocyto -c conda-forge numba pytables louvain

echo "scRNA_seurat_wrapper"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_seurat_wrapper/lib/R/library
./install/bin/conda create -n scRNA_seurat_wrapper -c r -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.6
./install/bin/conda install -n scRNA_seurat_wrapper -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr bioconductor-batchelor
#in R
/home/jichen/software/BETSY/install/envs/scRNA_seurat_wrapper/bin/R
install.packages("devtools")
library(devtools)
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
install.packages("dplyr")
devtools::install_github('satijalab/seurat-wrappers')



#scRNAcnv from BETSY pipeline
./install/bin/conda create -n scRNAcnv -c r -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.5.1 r-data.table r-matrix r-dplyr r-beeswarm r-knitr bioconductor-biomart
./install/envs/scRNAcnv/bin/R
install.packages("mclust")
#./install/bin/conda create -n scRNAcnv -c bioconda -c r r-base=3.5.1 bioconductor-biomart


#fishplot in a single env, 20190604
./install/bin/conda search r-base
./install/bin/conda create -n fishplot -c pkgs/main r-base=3.6.0
export R_LIBS=~/software/BETSY/install/envs/fishplot/lib/R/library
~/software/BETSY/install/envs/fishplot/bin/R
install.packages("devtools")
library(devtools)
install_github("chrisamiller/fishplot")

#pyclone in a single env, 20190604
./install/bin/conda create -n pyclone -c aroth85 pyclone
#install cluster 3.0 for pyclone
cd /home/jichen/software/BETSY/install/pkgs
wget http://bonsai.hgc.jp/~mdehoon/software/cluster/cluster-1.58.tar.gz
uz cluster-1.58.tar.gz
cd cluster-1.58
./configure --prefix=/home/jichen/software/BETSY/install/envs/pyclone/ --without-x
make
make install
#fix matplotlib DISPLAY errors
#https://github.com/aroth85/pyclone/issues/3
#change backend to "Agg" (Agg does not require X gui) in this file of your PYTHONPATH
vi /home/jichen/software/BETSY/install/envs/pyclone/lib/python2.7/site-packages/matplotlib/mpl-data/matplotlibrc
#backend      : qt5agg
backend      : Agg

conda activate pyclone
conda deactivate




#scRNA ssGSEA in a single env, 20190604
./install/bin/conda create -n ssGSEA -c bioconda bioconductor-gseabase
./install/bin/conda install -n ssGSEA -c bioconda -c anaconda python-gdbm bioconductor-gsva

#CNA in a single env, 20190605
./install/bin/conda create -n CNA -c conda-forge -c r -c bioconda sequenza-utils htslib
./install/bin/conda install -n CNA -c bioconda r-facets r-base=3.5.1

#Bulk RNA
./install/bin/conda create -n BulkRNAseq -c bioconda bioconductor-gseabase python-gdbm bioconductor-gsva bioconductor-edger bioconductor-s4vectors r-base=3.5.1
export R_LIBS=/home/jichen/software/BETSY/install/envs/BulkRNAseq/lib/R/library/
#DESeq2
/home/jichen/software/BETSY/install/envs/BulkRNAseq/bin/R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edger")
#font
cp ~/software/BETSY/genomicode/Verdana* /net/isi-dcnl/ifs/user_data/abild/jichen/software/Python_lib/lib/python2.7/site-packages/genomicode/
#povray
conda create -n BulkRNAseq_tools -c schrodinger povray=3.7.0.4
echo `pwd`/install/envs/BulkRNAseq_tools/bin/povray 
vi ~/.genomicoderc
#sva for combat
./install/bin/conda install -n BulkRNAseq -c bioconda bioconductor-sva

#WGS
 ./install/bin/conda create -n WGS -c bioconda qualimap


echo "Update BETSY for tools"
cd Update
python Update_genomicode.py --input /home/jichen/software/BETSY/install/bin/ --manual update_manual.list --output update20190523_2 > update20190523_2.log 2>&1 &
cp update20190523_2.genomicoderc ~/.genomicoderc
python Update_genomicode.py --input /home/jichen/software/BETSY/install/bin/ --manual update_manual.list --output update20190525_2 > update20190525_2.log 2>&1 &
cp update20190525_2.genomicoderc ~/.genomicoderc
python Update_genomicode.py --input /home/jichen/software/BETSY/install/bin/ --manual update_manual.list --output update20190528_2 > update20190528_2.log 2>&1 &
cp update20190528_2.genomicoderc ~/.genomicoderc


echo "install on Tikvah"
echo "snRNA in one place"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
export PATH=$PATH:/home/jichen/software/BETSY/install/envs/scRNA/bin/
#install seurat and singleR
./install/bin/conda create -n scRNA -c r -c anaconda -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.5.1 r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr r-devtools gxx_linux-64
./install/bin/conda install -n scRNA -c r r-devtools
./install/bin/conda install -n scRNA -c anaconda gxx_linux-64
./install/envs/scRNA/bin/R
library(devtools)
devtools::install_github('satijalab/seurat-wrappers')
devtools::install_github('dviraran/SingleR')

echo "infercnv"
./install/bin/conda create -n Infercnv -c pkgs/main -c r r-base=3.6.0
./install/bin/conda install -n Infercnv -c conda-forge jags
./install/bin/conda install -n Infercnv -c anaconda gxx_linux-64
export PKG_CONFIG_PATH=/home/jichen/software/BETSY/install/envs/Infercnv/lib/pkgconfig
export LD_RUN_PATH=/home/jichen/software/BETSY/install/envs/Infercnv/lib/JAGS/
export R_LIBS=/home/jichen/software/BETSY/install/envs/Infercnv/lib/R/library/
export PATH=$PATH:/home/jichen/software/BETSY/install/envs/Infercnv/bin/
./install/envs/Infercnv/bin/R
install.packages("rjags", configure.args="--enable-rpath")
library(rjags)
install.packages("BiocManager")
BiocManager::install("infercnv", dependencies=TRUE)
library(infercnv)

echo "SCG single cell DNA subclone"
#scg, https://bitbucket.org/aroth85/scg/wiki/Usage
./install/bin/conda create -n SCG -c anaconda pandas numpy scipy pyyaml python=2.7.6
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/SCG/lib/python2.7/site-packages/
python setup.py install --prefix /home/jichen/software/BETSY/install/envs/SCG/
#https://github.com/ucasdp/RobustClone R lib
./install/bin/conda install -n SCG -c r r=3.5.1
./install/bin/conda install -n SCG -c r -c conda-forge -c bioconda bioconductor-sva r-rann r-reshape r-shape r-gmodels r-igraph
./install/bin/conda install -n SCG -c r -c conda-forge -c bioconda bioconductor-genefilter r-vegan

echo "limma"
./install/bin/conda create -n limma -c r -c bioconda r=3.5.1 bioconductor-limma r-ggplot2 r-data.table r-dplyr bioconductor-biobase
./install/bin/conda install -n limma -c conda-forge r-ggpubr
export R_LIBS=/home/jichen/software/BETSY/install/envs/limma/lib/R/library/
./install/bin/conda install -n limma -c conda-forge r-ggalluvial

echo "RNAseq count normorlization and formating"
./install/bin/conda create -n count_transform -c r -c bioconda r=3.5.1 bioconductor-limma r-ggplot2 r-data.table r-dplyr r-tidyr r-reshape
export R_LIBS=/home/jichen/software/BETSY/install/envs/count_transform/lib/R/library/
#DESeq2
/home/jichen/software/BETSY/install/envs/count_transform/bin/R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("edger")


echo "scRNA_trajectory"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_trajectory/lib/R/library
./install/bin/conda create -n scRNA_trajectory -c r -c bioconda -c pkgs/r r-base=3.5.1
./install/bin/conda install -n scRNA_trajectory -c conda-forge udunits2 imagemagick
module load hdf5/1.10.3
module load singularity/3.4.2
export GITHUB_PAT='b2c67447730e038284f3a3e5fc85c16611d00fb5'
export SINGULARITY_CACHEDIR=/home/jichen/software/container/Singularity
/home/jichen/software/BETSY/install/envs/scRNA_trajectory/bin/R

echo "scRNA_enricher GSEA: clusterProfiler and DOSE"
./install/bin/conda create -n scRNA_enricher -c pkgs/main -c r r-base=3.6.0
/home/jichen/software/BETSY/install/envs/scRNA_enricher/bin/R
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
#read.gmt need gseabase
./install/bin/conda install -n scRNA_enricher -c bioconda bioconductor-gseabase bioconductor-gsva

echo "install dropbox"
./install/bin/conda create -n Dropbox -c anaconda dropbox

echo "install modin to speed up pandas"
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/Pandas_modin/lib/python3.6/site-packages/
./install/bin/conda create -n Pandas_modin python=3.6 pandas=0.23.4
conda activate Pandas_modin
python -m pip install modin[ray]
python -c "import modin.pandas as pd"

echo "install zinbwave for normalization"
./install/bin/conda create -n zinbwave -c conda-forge -c bioconda -c r r-ggplot2 bioconductor-zinbwave r-dplyr r-tidyr parallel bioconductor-scRNAseq r-tsne bioconductor-biomart r-umap
export R_LIBS=/home/jichen/software/BETSY/install/envs/zinbwave/lib/R/library

echo "breast cancer subtype: genefu"
export R_LIBS=/home/jichen/software/BETSY/install/envs/breast_cancer_subtype/lib/R/library
./install/bin/conda create -n breast_cancer_subtype -c pkgs/main -c r r-base=3.6.0
/home/jichen/software/BETSY/install/envs/breast_cancer_subtype/bin/R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("genefu")
install.packages("caret")
install.packages("xtable")

echo "install python3 and pyreadr for wring and reading Rdata and RDS"
#export PYTHONPATH=/home/jichen/software/BETSY/install/envs/python3/lib/python3.6/site-packages/
#./install/bin/conda create -n python3 -c conda-forge pyreadr
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/python_dask/lib/python3.8/site-packages
./install/bin/conda create -n python_dask -c conda-forge dask pyreadr


echo "install GSVA for ssGSEA using aritro's signature"
export R_LIBS=/home/jichen/software/BETSY/install/envs/GSVA_1.34/lib/R/library
./install/bin/conda create -n GSVA_1.34 -c bioconda bioconductor-gsva

echo "scRNA_MAST"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_MAST/lib/R/library
./install/bin/conda create -n scRNA_MAST -c r -c bioconda -c pkgs/r r-base=3.5.1
./install/bin/conda install -n scRNA_MAST -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr
/home/jichen/software/BETSY/install/envs/scRNA_MAST/bin/R
#R
install.packages("BiocManager")
BiocManager::install("MAST")
library(Seurat)
library(data.table)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(MAST)
install.packages("dplyr")

#R3.6
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA_MAST_R3.6/lib/R/library
./install/bin/conda create -n scRNA_MAST_R3.6 -c r -c bioconda -c pkgs/r r-base=3.6.0
./install/bin/conda install -n scRNA_MAST_R3.6 -c r -c bioconda -c conda-forge -c derickl r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr
/home/jichen/software/BETSY/install/envs/scRNA_MAST_R3.6/bin/R
install.packages("BiocManager")
BiocManager::install("MAST")

echo "DNA analysis"
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/lib/R/library:$R_LIBS
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/lib/python2.7/site-packages/$PYTHONPATH
export PERL5LIB=~/software/BETSY/install/envs/DNA_analysis_WES/lib/site_perl/5.26.2/:$PERL5LIB
export PATH=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/bin:$PATH

./install/bin/conda create -n DNA_analysis_WES -c r -c bioconda -c pkgs/r r-base python=2.7 vcf2maf
#in R
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
BiocManager::install("maftools")#
mkdir /home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib
cd /home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib
mkdir hg19
mkdir hg19/vep
../bin/vep_install.pl -a cf -s homo_sapiens -y GRCh37 -c /home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib/hg19/vep
#failed manal finish
export VEP_DATA=/home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib/hg19/vep/homo_sapiens/86_GRCh37
cd /home/jichen/software/BETSY/install/envs/DNA_analysis_WES/VEP_lib/hg19/vep/homo_sapiens/86_GRCh37
wget https://raw.githubusercontent.com/Ensembl/Bio-DB-HTS/master/scripts/convert_gz_2_bgz.sh
./convert_gz_2_bgz.sh $VEP_DATA/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz ~/software/BETSY/install/envs/DNA_analysis_WES/bin/bgzip


#This package installs only the variant effect predictor (VEP) library
#code. To install data libraries, you can use the 'vep_install.pl' command
#installed along with it. For example, to install the VEP library for human
#hg19/GRCh37 to a directory:

#vep_install.pl -a cf -s homo_sapiens -y GRCh37 -c /output/path/to/hg19/vep
#vep_convert_cache.pl -species homo_sapiens -version 86_GRCh37 -d /output/path/to/hg19/vep

#(note that vep_install.pl is renamed from INSTALL.pl
# and vep_convert_cache.pl from covert_cache.pl
# to avoid having generic script names in the PATH)

#The convert cache step is not required but improves lookup speeds during
#runs. See the VEP documentation for more details:

#http://useast.ensembl.org/info/docs/tools/vep/script/vep_cache.html

ech "install cloevlo"
./install/bin/conda create -n DNA_analysis_cloevol -c r -c bioconda -c pkgs/r r-base
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol/lib/R/library:$R_LIBS
/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol/bin/R
install.packages('devtools')
library(devtools)
install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
install_github("chrisamiller/fishplot")
./install/bin/conda install -n DNA_analysis_cloevol -c bioconda r-plotrix


echo "mutational pattern"
./install/bin/conda create -n DNA_mutational_pattern -c r -c bioconda -c pkgs/r r-base
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern/lib/R/library:$R_LIBS
/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern/bin/R
./install/bin/conda install -n DNA_mutational_pattern -c r r-nmf

echo "mutational pattern R3.6"
./install/bin/conda create -n DNA_mutational_pattern_R3.6 -c r -c bioconda -c pkgs/r r-base=3.6
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern_R3.6/lib/R/library
/home/jichen/software/BETSY/install/envs/DNA_mutational_pattern_R3.6/bin/R
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("MutationalPatterns")
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")
  BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
  install.packages("gridExtra")
  install.packages("ggplot2")
  install.packages("data.table")

echo "platypus for betsy germline variant calling"
./install/bin/conda create -n platypus-variant -c bioconda platypus-variant
#/net/isi-dcnl/ifs/user_data/abild/jichen/software/BETSY/install/envs/platypus-variant/share/platypus-variant-0.8.1.2-0/Platypus.py

echo "WES_cancer_cell_frequency"
./install/bin/conda create -n WES_cancer_cell_frequency -c r -c bioconda -c pkgs/r r-base
export R_LIBS=/home/jichen/software/BETSY/install/envs/WES_cancer_cell_frequency/lib/R/library:$R_LIBS
#./install/bin/conda install -n WES_cancer_cell_frequency -c anaconda gmp
./install/bin/conda install -n WES_cancer_cell_frequency -c conda-forge mpfr

echo "python2.7"
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/python2.7/lib/python2.7/site-packages/
./install/bin/conda create -n python2.7 -c bioconda pybedtools pandas numpy python=2.7 anaconda
./install/bin/conda install -n python2.7 -c conda-forge biopython
./install/bin/conda install -n python2.7 -c conda-forge pypdf2


echo "R tree: ggtree, ape"
export R_LIBS=/home/jichen/software/BETSY/install/envs/R_TREE/lib/R/library
./install/bin/conda create -n R_TREE -c bioconda -c conda-forge -c r r-base=3.5.1 bioconductor-ggtree r-ape r-phytools
#update the latest packakge to fix errors
#Error in UseMethod("full_join") : no applicable method for 'full_join' applied to an object of class "phylo"
/home/jichen/software/BETSY/install/envs/R_TREE/bin/R
BiocManager::install("treeio")
BiocManager::install("ggtree")
# R3.6
export R_LIBS=/home/jichen/software/BETSY/install/envs/R_TREE3.6/lib/R/library
./install/bin/conda create -n R_TREE3.6 -c bioconda -c conda-forge -c r r-base=3.6.1 bioconductor-ggtree bioconductor-treeio r-ape r-phytools r-tidyverse
## install cowplot to combine plots into one figure
./install/bin/conda install -n R_TREE3.6 -c conda-forge r-cowplot

echo "Infercnv_heatmap_R3.6"
export R_LIBS=/home/jichen/software/BETSY/install/envs/Infercnv_heatmap_R3.6/lib/R/library
./install/bin/conda create -n Infercnv_heatmap_R3.6 -c bioconda -c conda-forge -c r r-base=3.6.1 r-fastcluster r-data.table bioconductor-complexheatmap=2.0.0
./install/bin/conda install -n Infercnv_heatmap_R3.6 -c conda-forge r-dendextend r-tidyverse


echo "DENDRO_R3.6: scRNA phylogeny"
export R_LIBS=/home/jichen/software/BETSY/install/envs/DENDRO_R3.6/lib/R/library
./install/bin/conda create -n DENDRO_R3.6 -c r -c bioconda -c conda-forge r-base=3.6.1 r-data.table r-tidyverse
./install/bin/conda install -n DENDRO_R3.6 -c conda-forge r-haven
/home/jichen/software/BETSY/install/envs/DENDRO_R3.6/bin/R
BiocManager::install("Biobase")
devtools::install_github("zhouzilu/DENDRO")

echo "primer3"
./install/bin/conda create -n Primer3 -c bioconda primer3

echo "cardelino for clone ID"
export R_LIBS=/home/jichen/software/BETSY/install/envs/cardelino_R3.6/lib/R/library/
./install/bin/conda create -n cardelino_R3.6 -c r -c conda-forge -c bioconda r-devtools r-base=3.6.1 r-data.table r-tidyverse bioconductor-ggtree bioconductor-variantannotation r-biocmanager
./install/bin/conda install -n cardelino_R3.6 -c conda-forge r-foreign
/home/jichen/software/BETSY/install/envs/cardelino_R3.6/bin/R
#devtools::install_github("PMBio/cardelino", build_vignettes = TRUE)

echo "monocle3"
export R_LIBS=/home/jichen/software/BETSY/install/envs/monocle3_R3.6/lib/R/library/
./install/bin/conda create -n monocle3_R3.6 -c r -c anaconda -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.6.1 r-data.table r-matrix r-dplyr r-seurat r-beeswarm umap-learn=0.3.8 r-knitr r-devtools gxx_linux-64 r-biocmanager
/home/jichen/software/BETSY/install/envs/monocle3_R3.6/bin/R
# monocle2 installled, v3 failed
BiocManager::install(c("monocle"))
devtools::install_github(repo = 'satijalab/seurat', ref = 'develop')
./install/bin/conda install -n monocle3_R3.6 -c bioconda -c conda-forge r-glue=1.3.2

echo "install monocle3 beta, https://cole-trapnell-lab.github.io/monocle3/docs/installation/"
# monocle2 can not handle too many cells
#> monocle_cds <- orderCells(monocle_cds)
#Error in graph.adjacency.dense(adjmatrix, mode = mode, weighted = weighted,  : 
#  long vectors not supported yet: ../../src/include/Rinlinedfuns.h:522
export R_LIBS=/home/jichen/software/BETSY/install/envs/monocle3beta_R3.6/lib/R/library/
./install/bin/conda create -n monocle3beta_R3.6 -c r -c anaconda -c bioconda -c pkgs/r -c conda-forge -c derickl r-base=3.6.1
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge udunits2
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge r-sf
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge r-spdep
/home/jichen/software/BETSY/install/envs/monocle3beta_R3.6/bin/R
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/leidenbase')
devtools::install_github('cole-trapnell-lab/monocle3')
#Fix error
#monocle_cds <- cluster_cells(monocle_cds, reduction_method="UMAP")
#Error in dyn.load(file, DLLpath = DLLpath, ...) : 
#  unable to load shared object '/net/isi-dcnl/ifs/user_data/abild/jichen/software/BETSY/install/envs/monocle3beta_R3.6/lib/R/library/leidenbase/libs/leidenbase.so':
#  libRlapack.so: cannot open shared object file: No such file or directory
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge lapack
devtools::install_github('cole-trapnell-lab/leidenbase')
##install seuratv3
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge r-fitdistrplus
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge r-leiden
./install/bin/conda install -n monocle3beta_R3.6 -c conda-forge r-seurat

echo "cell cycle phase: reCAT"
export R_LIBS=/home/jichen/software/BETSY/install/envs/reCAT_R3.6/lib/R/library/
 1213  ./install/bin/conda create -n reCAT_R3.6 -c r -c conda-forge r-tsp r-ggplot2 r-cluster r-mclust r-doparallel

echo "pheatmap"
export R_LIBS=/home/jichen/software/BETSY/install/envs/pheatmap/lib/R/library/
./install/bin/conda create -n pheatmap -c r -c bioconda -c conda-forge bioconductor-limma r-ggplot2 r-ggpubr r-base=3.6.1 r-data.table r-tidyverse

echo "phyloWGS"
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/phyloWGS/lib/python2.7/site-packages/
export R_LIBS=/home/jichen/software/BETSY/install/envs/phyloWGS/lib/R/library/
./install/bin/conda create -n phyloWGS -c bioconda -c etetoolkit -c conda-forge scipy pandas numpy python=2.7 ete2
./install/bin/conda install -n phyloWGS -c bioconda phylowgs
./install/bin/conda install -n phyloWGS -c conda-forge gsl
./install/bin/conda install -n phyloWGS -c conda-forge pandas
export LD_LIBRARY_PATH=/home/jichen/software/BETSY/install/envs/phyloWGS/lib/:$LD_LIBRARY_PATH
cd /home/jichen/software/Variants_calling/Subclone_tools/PhyloWGS/phylowgs/ 
g++ -o mh.o -O3 mh.cpp  util.cpp `gsl-config --cflags --libs`
cp mh.o /home/jichen/software/BETSY/install/envs/phyloWGS/bin/
#cp other python libraries too (util2 et al).
./install/bin/conda install -n phyloWGS -c conda-forge r-devtools
library(devtools)
install_github('MathOnco/EvoFreq')

echo "R_combine_pdf_figures"
export R_LIBS=/home/jichen/software/BETSY/install/envs/R_combine_pdf_figures/lib/R/library/
./install/bin/conda create -n R_combine_pdf_figures -c conda-forge r-pdftools r-magick
./install/bin/conda install -n R_combine_pdf_figures -c conda-forge -c r r-ggplot2 r-cowplot r-tidyverse

ech "install cloevlo"
./install/bin/conda create -n DNA_analysis_cloevol_R -c r -c bioconda -c pkgs/r r-base
export R_LIBS=/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/lib/R/library
/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/R
install.packages('devtools')
library(devtools)
install_github('hdng/clonevol')
install.packages('gridBase')
install.packages('gridExtra')
install.packages('ggplot2')
install.packages('igraph')
install.packages('packcircles')
install_github('hdng/trees')
./install/bin/conda install -n DNA_analysis_cloevol_R -c conda-forge r-foreign
/home/jichen/software/BETSY/install/envs/DNA_analysis_cloevol_R/bin/R
install_github("chrisamiller/fishplot")
./install/bin/conda install -n DNA_analysis_cloevol_R -c bioconda r-plotrix r-data.table
./install/bin/conda install -n DNA_analysis_cloevol_R -c r r-reshape2
./install/bin/conda install -n DNA_analysis_cloevol_R -c r r-dplyr
./install/bin/conda install -n DNA_analysis_cloevol_R -c conda-forge r-rlang

echo "MEDALT_sc_tree"
export R_LIBS=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/lib/R/library
export PATH=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/:$PATH
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/lib/python2.7/site-packages/
./install/bin/conda create -n MEDALT_sc_tree -c r -c bioconda -c pkgs/r -c conda-forge python-igraph r-base bioconductor-helloranges
./install/bin/conda install -n MEDALT_sc_tree -c conda-forge r-devtools
./install/bin/conda install -n MEDALT_sc_tree -c conda-forge r-rlang r-mass
./install/bin/conda install -n MEDALT_sc_tree -c conda-forge r-circlize
/home/jichen/software/BETSY/install/envs/MEDALT_sc_tree/bin/R

#echo "MEDALT_sc_tree_R3.5"
#./install/bin/conda create -n MEDALT_sc_tree_R3.5 -c r -c bioconda -c pkgs/r -c conda-forge python-igraph r-base=3.5 bioconductor-helloranges r-devtools r-rlang r-mass
#./install/bin/conda install -n MEDALT_sc_tree_R3.5 -c r r-hmisc


echo "ABSOLUTE"
export R_LIBS=/home/jichen/software/BETSY/install/envs/ABSOLUTE/lib/R/library
./install/bin/conda create -n ABSOLUTE -c astro r devtools
cd ~/software/Variants_calling/Subclone_tools/
#install numDeriv
#/home/jichen/software/BETSY/install/envs/ABSOLUTE/bin/R
#install.packages("numDeriv")
/home/jichen/software/BETSY/install/envs/ABSOLUTE/bin/R CMD INSTALL ABSOLUTE_1.0.6.tar.gz
#install doabsolute
#/home/jichen/software/BETSY/install/envs/ABSOLUTE/bin/R
#devtools::install_github("ShixiangWang/DoAbsolute")
./install/bin/conda install -n ABSOLUTE -c r r-tidyverse

echo "graphlan: high-quality circular representations of taxonomic and phylogenetic trees"
./install/bin/conda create -n graphlan -c bioconda graphlan
export PATH=/home/jichen/software/BETSY/install/envs/graphlan/bin/:$PATH
export PYTHONPATH=/home/jichen/software/BETSY/install/envs/graphlan/lib/python2.7/site-packages/

