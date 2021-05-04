#!/usr/bin/env R

# Test IntEREst tool using provided vignette (https://www.bioconductor.org/packages/release/bioc/vignettes/IntEREst/inst/doc/IntEREst.html)
# 
# Open R after activating the environment using:
# 
# > conda activate r403_interest
# 
# 

library(IntEREst)

#---------------------
# prep annotation data
#---------------------

tmpGen<-u12[u12[,"gene_name"]=="RHBDD3",] # data for single gene
dim(tmpGen)
# [1] 15 16
class(tmpGen)
# [1] "data.frame"
colnames(tmpGen)
# [1] "id"           "int_ex_id"    "chr"          "begin"        "end"
# [6] "strand"       "int_ex"       "trans_type"   "ens_gene_id"  "ens_trans_id"
# [11] "int_ex_num"   "gene_name"    "trans_name"   "overlap_no"   "int_type"
# [16] "int_subtype"

tmpEx<-tmpGen[tmpGen[,"int_ex"]=="exon",] # get gene exons
dim(tmpEx)
# [1]  8 16
class(tmpEx)
# [1] "data.frame"

# build gff3 file
exonDat<- cbind(tmpEx[,3], ".",
    tmpEx[,c(7,4,5)], ".", tmpEx[,6], ".",paste("ID=exon",
    tmpEx[,11], "; Parent=ENST00000413811", sep="") )
class(exonDat)
# [1] "data.frame"
dim(exonDat)
# [1] 8 9

trDat <- c(tmpEx[1,3], ".", "mRNA", as.numeric(min(tmpEx[,4])),
    as.numeric(max(tmpEx[,5])), ".", tmpEx[1,6], ".",
    "ID=ENST00000413811")
class(trDat)
# [1] "character"
trDat
# [1] "chr22"              "."                  "mRNA"
# [4] "29655841"           "29664198"           "."
# [7] "-"                  "."                  "ID=ENST00000413811"

outDir<- file.path(tempdir(),"tmpFolder"); dir.create(outDir)
outDir<- normalizePath(outDir); outDir
# "/private/var/folders/jf/8cq6vdx54rb7hlwv4d7zwmmr5s1kwv/T/RtmpSHnCXA/tmpFolder"

gff3File<-paste(outDir, "gffFile.gff", sep="/")
cat("##gff-version 3\n",file=gff3File, append=FALSE)
cat(paste(paste(trDat, collapse="\t"),"\n", sep=""), file=gff3File, append=TRUE)
write.table(exonDat, gff3File, row.names=FALSE, col.names=FALSE, 
    sep='\t', quote=FALSE, append=TRUE)

# get u12 introns from data
u12Int<-u12[u12$int_ex=="intron"&u12$int_type=="U12",]
class(u12Int)
# [1] "data.frame"
dim(u12Int)
# [1] 555  16

# build the reference
testRef <- referencePrepare(sourceBuild="file", filePath=gff3File, 
    u12IntronsChr=u12Int[,"chr"], u12IntronsBeg=u12Int[,"begin"], 
    u12IntronsEnd=u12Int[,"end"], collapseExons=TRUE, fileFormat="gff3", 
    annotateGeneIds=FALSE)
class(testRef)
# [1] "data.frame"
dim(testRef)
# [1] 15  9
colnames(testRef)
# [1] "chr"                      "begin"
# [3] "end"                      "strand"
# [5] "int_ex"                   "int_ex_num"
# [7] "collapsed_transcripts_id" "collapsed_transcripts"
# [9] "int_type"

#----------------------
# get human genome hg19
#----------------------
# install genome bioc library
# note: large file download (677.3 MB)
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# import the genome
BSgenome.Hsapiens.UCSC.hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19



