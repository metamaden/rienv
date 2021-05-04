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
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.UCSC.hg19")

# import the genome
BSgenome.Hsapiens.UCSC.hg19 <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

# subset genes
ind<- u12$gene_name %in% c("RHBDD2", "YBX2")

# annotate u12 data

pwm.u12 <- list(pwmU12db[[1]][,11:17],pwmU12db[[2]],
        pwmU12db[[3]][,38:40],pwmU12db[[4]][,11:17],
        pwmU12db[[5]][,38:40])
# pwm.u12 description:
# A list containing position weight matrices of (in order): Donor site, branch 
# point, and acceptor site of U12-type introns, and donor site and acceptor 
# site of U2-type introns. If not provided, the information related to pwmU12db 
# data is used.

pwm.ssindex <- list(indexDonU12=1, indexBpU12=1, indexAccU12=3,
        indexDonU2=1, indexAccU2=3)
# pwm.ssindex description:
# A list (or vector) that contains the column number in each element of 
# pwmU12U2 that represents the 5' or 3' Splice Site; The order should be 
# equivalent to the pwmU12U2. If not provided the information from pwmU12db 
# data is used, i.e. pwmSsIndex=list(indexDonU12=1, indexBpU12=1, 
# indexAccU12=3, indexDonU2=1, indexAccU2=3)

min.match.score <- c(rep(paste(80,"%",sep=""),2), "40%", paste(80,"%",sep=""), "40%")
# min.match.score description:
# Min percentage match score, when scoring matching of a sequence to pwm. 
# Different score thresholds could also be defined for the various sites 
# (U12/U2 donors, the U12 branch point and U12/U2 acceptors); A vector with 
# 5 elements can be assigned which each shows the match score to use for each 
# PWM in pwmU12U2.

annoU12 <- annotateU12(pwmU12U2 = pwm.u12, pwmSsIndex = pwm.ssindex,
    referenceChr = u12[ind,'chr'], referenceBegin = u12[ind,'begin'],
    referenceEnd = u12[ind,'end'], referenceIntronExon = u12[ind,"int_ex"],
    intronExon = "intron", matchWindowRelativeUpstreamPos = c(NA,-29,NA,NA,NA),
    matchWindowRelativeDownstreamPos=c(NA,-9,NA,NA,NA),
    minMatchScore = min.match.score, refGenome = BSgenome.Hsapiens.UCSC.hg19, 
    setNaAs = "U2", annotateU12Subtype = TRUE)
# annotateU12 description:
# Receives coordinates, a reference genome and PWMs of splice site of U12 and 
# U2 type introns, and returns a data.frame with 2 columns. The first column 
# shows wheather the corresponding sequences matches U12, U2 or both (U12/U2) 
# consensus sequences (based on their score when fitting the PWMs). The second 
# column shows whether the match is on positive strand or negative when fitting 
# the PWMs to the sequences

# returns: 
## 
## U12  U2 
##   2  10
# equivalent to: table(annoU12[,1])

head(annoU12)
#  int_type strandMatch u12_subtype
#1     <NA>        <NA>        <NA>
#2       U2        <NA>        <NA>
#3     <NA>        <NA>        <NA>
#4       U2        <NA>        <NA>
#5     <NA>        <NA>        <NA>
#6      U12           +       GT-AG

#------------------------------
# intron retention calculations
#------------------------------
# manage output path
outDir<- file.path(tempdir(),"interestFolder"); dir.create(outDir); 
outDir <- normalizePath(outDir)

# load bam file
bamF <- system.file("extdata", "small_test_SRR1691637_ZRSR2Mut_RHBDD3.bam",
    package="IntEREst", mustWork=TRUE)

# choose reference for gene
ref<-u12[u12[,"gene_name"]=="RHBDD3",]

#------------
# IR analysis
#------------
# Reads mapping to inner introns are considered, hence junctionReadsOnly is FALSE
testInterest <- interest(bamFileYieldSize = 10000, junctionReadsOnly = FALSE,
    bamFile = bamF, isPaired = TRUE, isPairedDuplicate = FALSE,
    isSingleReadDuplicate = NA, reference = ref, 
    referenceGeneNames = ref[,"ens_gene_id"], referenceIntronExon = ref[,"int_ex"],
    repeatsTableToFilter = c(), outFile = paste(outDir, "intRetRes.tsv", sep="/"),
    logFile = paste(outDir, "log.txt", sep="/"), method = "IntRet", clusterNo = 1, 
    returnObj = FALSE, scaleLength = TRUE, scaleFragment = TRUE)

# interest() description:
# A read summarization function that countsns all the reads mapping to the 
# introns/exons based on the users detailed parameter settings. The process 
# can be run in parallel on multiple computing cores to improve it performance.

# returns:
# InERESt: Running bamPreprocess. Detailed log info are written in:  
# /private/var/folders/jf/8cq6vdx54rb7hlwv4d7zwmmr5s1kwv/T/RtmpSHnCXA/interestFolder/log.txt
# Log info: Running interest in Parallel mode.
# InERESt: Running interestAnalyse.
# InERESt:interestAnalyse: Begins ...
# Analyzing paired reads, ID:  SRR1691633.22391522
# Analyzing single reads, ID:  SRR1691633.10002641
# InERESt:interestAnalyse: Read counting ends. Running time:  0.97349  secs
# InERESt: Running interestSummarise.
# IntERESt:interestSummarise: Begins ...
# IntERESt:interestSummarise: Normalizing intron retention read levels.
# InERESt: run ends. Full running time:  0.980072  secs

testIntRetObj<- readInterestResults(resultFiles = paste(outDir, "intRetRes.tsv", sep="/"), 
    sampleNames="small_test_SRR1691637_ZRSR2Mut_RHBDD3", 
    sampleAnnotation=data.frame(type="ZRSR2mut", test_ctrl="test"), 
    commonColumns=1:ncol(ref), freqCol=ncol(ref)+1, 
    scaledRetentionCol=ncol(ref)+2, scaleLength=TRUE, scaleFragment=TRUE, 
    reScale=TRUE, geneIdCol="ens_gene_id")

# readInterestResults() description:
# Reads one or multiple text file results generated by the interest or 
# interest.sequential functions and builds an object of 
# SummarizedExperiment-class class.

testIntRetObj
# class: SummarizedExperiment
# dim: 15 1
# metadata(2): scaleFragment scaleLength
# assays(2): counts scaledRetention
# rownames: NULL
# rowData names(16): id int_ex_id ... int_type int_subtype
# colnames(1): small_test_SRR1691637_ZRSR2Mut_RHBDD3
# colData names(3): resultFiles type test_ctrl

# returns:
# Reading file 1 / 1 : /private/var/folders/jf/8cq6vdx54rb7hlwv4d7zwmmr5s1kwv/
# T/RtmpSHnCXA/interestFolder/intRetRes.tsv

# Reads mapping to inner introns are considered, hence junctionReadsOnly is FALSE
testInterest <- interest(bamFileYieldSize=10000, junctionReadsOnly=FALSE,
    bamFile=bamF, isPaired=TRUE, isPairedDuplicate=FALSE, isSingleReadDuplicate=NA,
    reference=ref, referenceGeneNames=ref[,"ens_gene_id"], 
    referenceIntronExon=ref[,"int_ex"], repeatsTableToFilter=c(),
    outFile=paste(outDir, "intSpanRes.tsv", sep="/"), 
    logFile=paste(outDir, "log.txt", sep="/"), method="IntSpan", clusterNo=1, 
    returnObj=FALSE, scaleLength= TRUE, scaleFragment= TRUE)

# returns:
# InERESt: Running bamPreprocess. Detailed log info are written in:  /private/var/folders/jf/8cq6vdx54rb7hlwv4d7zwmmr5s1kwv/T/RtmpSHnCXA/interestFolder/log.txt
# Log info: Running interest in Parallel mode.
# InERESt: Running interestAnalyse.
# InERESt:interestAnalyse: Begins ...
# Analyzing paired reads, ID:  SRR1691633.22391522
# Analyzing single reads, ID:  SRR1691633.10002641
# InERESt:interestAnalyse: Read counting ends. Running time:  0.6759219  secs
# InERESt: Running interestSummarise.
# IntERESt:interestSummarise: Begins ...
# IntERESt:interestSummarise: Normalizing intron retention read levels.
# InERESt: run ends. Full running time:  0.681663  secs

testIntSpanObj<- readInterestResults(resultFiles= paste(outDir, "intSpanRes.tsv", sep="/"), 
    sampleNames="small_test_SRR1691637_ZRSR2Mut_RHBDD3", 
    sampleAnnotation=data.frame(type="ZRSR2mut", test_ctrl="test"), 
    commonColumns=1:ncol(ref), freqCol=ncol(ref)+1, scaledRetentionCol=ncol(ref)+2, 
    scaleLength=TRUE, scaleFragment=TRUE, reScale=TRUE, geneIdCol="ens_gene_id")

# returns:
# Reading file 1 / 1 : /private/var/folders/jf/8cq6vdx54rb7hlwv4d7zwmmr5s1kwv/T
# /RtmpSHnCXA/interestFolder/intSpanRes.tsv

testIntSpanObj
# class: SummarizedExperiment
# dim: 15 1
# metadata(2): scaleFragment scaleLength
# assays(2): counts scaledRetention
# rownames: NULL
# rowData names(16): id int_ex_id ... int_type int_subtype
# colnames(1): small_test_SRR1691637_ZRSR2Mut_RHBDD3
# colData names(3): resultFiles type test_ctrl

# get exon data
testExExObj <- readInterestResults(resultFiles= paste(outDir, "exExRes.tsv", sep="/"), 
    sampleNames="small_test_SRR1691637_ZRSR2Mut_RHBDD3", 
    sampleAnnotation=data.frame(type="ZRSR2mut", test_ctrl="test"), 
    commonColumns=1:ncol(ref), freqCol=ncol(ref)+1, scaledRetentionCol=ncol(ref)+2, 
    scaleLength=TRUE, scaleFragment=TRUE, reScale=TRUE, geneIdCol="ens_gene_id")

testExExObj

#---------------
# compare counts
#---------------
head(counts(testIntRetObj))
# small_test_SRR1691637_ZRSR2Mut_RHBDD3
# [1,]                                     0
# [2,]                                    11
# [3,]                                     0
# [4,]                                    11
# [5,]                                     0
# [6,]                                   180

head(counts(testIntSpanObj))
head(counts(testExExObj))