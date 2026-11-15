# Optional local fixture using the same featureCounts engine via Rsubread.
args <- commandArgs(trailingOnly=TRUE)
directory <- args[1]; dir.create(directory,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(Rsubread))
annotation <- data.frame(GeneID=c("A","A","B","C"),Chr="chr1",Start=c(100,200,250,100),End=c(150,300,350,200),Strand=c("+","+","+","-"))
records <- c("@HD\tVN:1.6\tSO:coordinate","@SQ\tSN:chr1\tLN:1000")
for (i in 1:6) {
  position <- c(120,120,140,280,320,500)[i]
  cigar <- if(i==3) "11M49N20M" else "20M"
  length <- if(i==3) 31 else 20
  records <- c(records,paste(paste0("r",i),if(i==2) 0 else 16,"chr1",position,255,cigar,"*",0,0,
    paste(rep("A",length),collapse=""),paste(rep("I",length),collapse=""),"NH:i:1",sep="\t"))
}
sam <- file.path(directory,"fixture.sam"); writeLines(records,sam)
result <- featureCounts(sam,annot.ext=annotation,isGTFAnnotationFile=FALSE,strandSpecific=2,
                       allowMultiOverlap=FALSE,countMultiMappingReads=FALSE,primaryOnly=TRUE,nthreads=1)
stopifnot(identical(as.integer(result$counts[c("A","B","C"),1]),c(2L,1L,1L)))
stopifnot(result$stat[result$stat$Status=="Unassigned_Ambiguity",2]==1,
          result$stat[result$stat$Status=="Unassigned_NoFeatures",2]==1)
cat("PASS: reverse-stranded SAF; same-gene exons count once; opposite strand is separate; multi-gene ambiguity excluded.\n")
