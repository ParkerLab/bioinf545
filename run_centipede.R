#!/usr/bin/Rscript
# ---------------------------------------------------------------------------


library(CENTIPEDE)

args<-commandArgs(TRUE)

# args1 = Tab-delimited gzipped cuts file (Required)
# args2 = Tab-delimited gzipped motif file with score in the 5th column (Required)
# args3 = name for tab-delimited gzipped output file that contains CENTIPEDE posterior probability appended to motifs file (Required)
# args4 = column that contain motif scores. Defaults to 5
# args5 = number of lines to read in. Defaults to -1

scores <- ifelse(!is.na(args[4]), as.numeric(args[4]),  5)
nlines <- ifelse(!is.na(args[5]), as.numeric(args[5]), -1)


###################### Read Data & Manipulate ################################

cuts<-read.table(file=args[1], colClasses = 'integer', nrows = nlines)

anno<-read.table(file=args[2], nrows = nlines)

centFit <-fitCentipede(Xlist = list(DNase=as.matrix(cuts)), Y=cbind(rep(1, dim(anno)[1]), anno[,scores]), DampLambda=0.01,DampNegBin=0.001)

output<-cbind(anno, centFit$PostPr)

write.table(output, gzfile(args[3]), quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

