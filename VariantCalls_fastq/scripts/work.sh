
echo "WES fastq file list"
#FELINE 13,15,19,20,21,22,25,28,29,31,34,35,36,40,41,43
ls /home/jichen/abild/U54_data_folder/COH068/17639R/Fastq/*.fastq.gz >> file_fq.list
cp /home/jichen/abild/U54_data_folder/COH068/feline_47.sample.txt ./
#FELINE 14,23,24,26,27,37,39
ls /home/jichen/abild/U54_data_folder/COH077/FT-SH6178/*/*.fastq.gz >> file_fq.list
xls2txt /home/jichen/abild/U54_data_folder/COH077/feline_7patients.xlsx > feline_7patients.txt
#FELINE 44,45,46
ls /home/jichen/abild/U54_data_folder/COH074/FT-SH6014/*/*.gz >> file_fq.list
#There are some additional WES of FELINE 25, 31, 43 in COH079. Add these sequences make the FACETS/sequenza calls weired. I did not include these in my analysis.

echo "merge WES fastq info"
export R_LIBS=/home/jichen/software/BETSY/install/envs/scRNA/lib/R/library
/home/jichen/software/BETSY/install/envs/scRNA/bin/Rscript FELINE_WES_fastq.R
grep -v "Sample" FELINE_WES_fastq.list | sort -k1,1 -k4,4 -k2,2 >> FELINE_WES_fastq.sorted.list

