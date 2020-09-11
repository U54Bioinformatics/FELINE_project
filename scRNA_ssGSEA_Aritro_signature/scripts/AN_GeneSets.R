load("AN_GeneSets.RData")
lapply(AN_GeneSets, write, "AN_GeneSets.gmt", append=TRUE, ncolumns=1000, sep="\t")

