# Offline database fixture; actual hypergeometric/BH/fgsea and plot code runs.
args <- commandArgs(trailingOnly=TRUE)
root <- args[1]; dir.create(root,recursive=TRUE,showWarnings=FALSE)
script <- args[2]
genes <- data.frame(gene_id=paste0("g",1:600),background=1,foreground=as.integer(1:600 <= 120),
  direction=ifelse(1:600<=60,"any_dge;upregulated",ifelse(1:600<=120,"any_dge;downregulated","")),
  rank_score=c(seq(12,5,length.out=60),seq(-12,-5,length.out=60),seq(-4,4,length.out=480)))
input <- file.path(root,"prepared.tsv")
write.table(genes,input,sep="\t",quote=FALSE,row.names=FALSE)
arguments <- c("--input",input,"--outdir",root,"--analysis-type","dge","--species","mouse",
               "--genome","GRCm39","--contrast-id","B_vs_A","--numerator","B","--denominator","A",
               "--go","false","--reactome","false","--kegg","false","--hallmarks","true","--rich-plots","false")
environment <- new.env(parent=globalenv())
environment$commandArgs <- function(trailingOnly=FALSE) arguments
for (expression in parse(script)) {
  eval(expression,envir=environment)
  if(is.call(expression) && identical(expression[[1]],as.name("<-")) && identical(expression[[2]],as.name("msig_call"))) {
    environment$msig_call <- function(collection,subcollection=NULL,database) {
      data.frame(database=database,term_id=rep(c("positive","negative","neutral"),each=30),
        term_name=rep(c("Positive fixture","Negative fixture","Neutral fixture"),each=30),
        gene_id=paste0("g",c(1:30,61:90,301:330)),db_version="offline_fixture_v1")
    }
  }
}
gsea <- read.delim(file.path(root,"gsea.tsv"))
ora <- read.delim(file.path(root,"ora.tsv"))
plots <- read.delim(file.path(root,"plot_index.tsv"))
stopifnot(gsea$NES[gsea$term_id=="positive"]>0,gsea$NES[gsea$term_id=="negative"]<0,
          all(c("upregulated","downregulated","any_dge") %in% ora$query),
          all(c("positive","negative") %in% plots$query[plots$method=="GSEA"]),
          all(file.exists(plots$png)),all(file.exists(plots$pdf)),file.exists(file.path(root,"mapping_genes.tsv")))
cat("PASS: real ORA/BH/fgsea, separate methods and directions, signed NES, per-gene mapping and plots; offline fixture database.\n")
