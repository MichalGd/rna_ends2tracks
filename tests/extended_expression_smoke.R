# Real DESeq2/edgeR fixture; run with an isolated output directory argument.
args <- commandArgs(trailingOnly=TRUE)
root <- normalizePath(args[1],winslash="/",mustWork=FALSE)
script <- normalizePath(args[2],winslash="/",mustWork=TRUE)
dir.create(root,recursive=TRUE,showWarnings=FALSE)
set.seed(13)
depth <- c(.5,1,2,1,.6,1.9)
mu <- rep(exp(runif(600,log(30),log(1200))),6)*rep(depth,each=600)
matrix <- matrix(rnbinom(3600,mu=mu,size=25),nrow=600)
matrix[1:60,4:6] <- matrix[1:60,4:6]*6L
matrix[61:120,1:3] <- matrix[61:120,1:3]*6L
matrix <- rbind(matrix,one_condition=c(0,0,0,20,20,40),single_sample=c(100,0,0,0,0,0),too_low=rep(1,6))
rownames(matrix) <- c(paste0("g",1:600),"one_condition","single_sample","too_low")
colnames(matrix) <- paste0("s",1:6)
write_tsv <- function(x,path) write.table(x,path,sep="\t",quote=FALSE,row.names=FALSE)
write_tsv(data.frame(gene_id=rownames(matrix),matrix,check.names=FALSE),file.path(root,"counts.tsv"))
write_tsv(data.frame(sample_id=colnames(matrix),condition=rep(c("A","B"),each=3)),file.path(root,"samples.tsv"))
write_tsv(data.frame(contrast_id="B_vs_A",factor="condition",numerator="B",denominator="A",
  design_mode="unpaired",resolved_design="~ condition",paired=FALSE,n_pairs=0,design_status="unpaired"),file.path(root,"contrasts.tsv"))
rscript <- file.path(R.home("bin"),"Rscript")
run <- function(mode,extra=character()) {
  command <- c("--vanilla",script,"--counts",file.path(root,"counts.tsv"),"--samples",file.path(root,"samples.tsv"),
    "--contrasts",file.path(root,"contrasts.tsv"),"--outdir",file.path(root,"results"),"--design","~ condition",
    "--mode",mode,"--contrast-id","B_vs_A", "--index-file",file.path(root,"index.tsv"),extra)
  status <- system2(rscript,vapply(command,shQuote,character(1)),stdout=file.path(root,paste0(mode,".log")),stderr="")
  if(status != 0) stop(paste("Fixture failed; inspect",root,mode))
}
run("qc")
run("contrast")
directory <- file.path(root,"results","B_vs_A")
read <- function(name) read.delim(file.path(directory,name),check.names=FALSE)
strict <- read("stringent/deseq2.tsv"); relaxed <- read("relaxed_fdr/deseq2.tsv")
stopifnot(identical(strict$gene_id,relaxed$gene_id),identical(strict$valid_test,relaxed$valid_test))
stopifnot(sum(strict$padj <= .05,na.rm=TRUE) > 0, sum(strict$padj <= .05,na.rm=TRUE) <= sum(relaxed$padj <= .05,na.rm=TRUE))
audit <- read("eligibility.tsv")
stopifnot(audit$eligible[audit$gene_id == "one_condition"], !audit$eligible[audit$gene_id == "single_sample"], !audit$eligible[audit$gene_id == "too_low"])
ranking <- read("ranking.tsv")
stopifnot(identical(ranking$stat,relaxed$stat))
model <- tools::md5sum(file.path(directory,"model.rds")); ranking_hash <- tools::md5sum(file.path(directory,"ranking.tsv"))
run("contrast",c("--stringent","false"))
index <- read.delim(file.path(root,"index.tsv"))
stopifnot(!"stringent" %in% index$analysis_mode, identical(model,tools::md5sum(file.path(directory,"model.rds"))),
          identical(ranking_hash,tools::md5sum(file.path(directory,"ranking.tsv"))))
cat("PASS: real QC/DESeq2/edgeR, unequal depths, one-condition support, common universe, formal thresholds, shared ranks, cached fit, switch off.\n")
cat("Eligible:",sum(audit$eligible),"Stringent:",sum(strict$padj<=.05,na.rm=TRUE),"Relaxed:",sum(relaxed$padj<=.05,na.rm=TRUE),"\n")
