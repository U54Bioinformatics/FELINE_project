outfile=FELINE.pyclone_mutation_num_sequenza.20_5_0.05_1_cnvall.txt
if [ ! -e $outfile ]; then
   echo "sum sequenza mutation number"
   cat FEL0*_pyclone_mutation_sequenza/num_mutations.txt | head -n 1 > $outfile
   cat FEL0*_pyclone_mutation_sequenza/num_mutations.txt | grep -v "Orig"  >> $outfile
fi

outfile=FELINE.pyclone_mutation_num_FACETS.20_5_0.05_1_cnvall.txt
if [ ! -e $outfile ]; then
   echo "sum FACETS mutation number"
   cat FEL0*_pyclone_mutation_FACETS/num_mutations.txt | head -n 1 > $outfile
   cat FEL0*_pyclone_mutation_FACETS/num_mutations.txt | grep -v "Orig"  >> $outfile
fi
