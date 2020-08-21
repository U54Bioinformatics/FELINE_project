args <- commandArgs(trailingOnly = TRUE)
infile=args[1]
outfile=args[2]

x = read.table(infile,sep="\t")
write.table(x, outfile, quote=F, row.names=F, sep="\t")

