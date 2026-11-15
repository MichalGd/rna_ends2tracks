args <- commandArgs(trailingOnly=TRUE)
arg <- function(key, default=NULL) {
  i <- match(key, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing value for", key))
  args[i+1L]
}
suppressPackageStartupMessages({library(DESeq2); library(edgeR); library(ggplot2)})
sample_distance_plot <- function(distance_matrix) {
  # Cluster the original VST distances, not distances between distance profiles.
  tree <- hclust(as.dist(distance_matrix), method="complete")
  sample_order <- rownames(distance_matrix)[tree$order]
  long <- data.frame(
    sample_x=factor(rep(colnames(distance_matrix), each=nrow(distance_matrix)), levels=sample_order),
    sample_y=factor(rep(rownames(distance_matrix), times=ncol(distance_matrix)), levels=rev(sample_order)),
    distance=as.vector(distance_matrix)
  )
  ggplot(long, aes(sample_x, sample_y, fill=distance)) + geom_tile() +
    scale_fill_gradient(low="white", high="#2166AC") + coord_fixed() +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) + theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1), panel.grid=element_blank()) +
    labs(x=NULL, y=NULL, fill="VST distance", title="VST sample distances",
         subtitle="Euclidean distance; complete-linkage hierarchical clustering")
}
required <- c("--counts", "--samples", "--contrasts", "--outdir", "--design")
for (key in required) if (is.null(arg(key))) stop(paste("Required argument:", key))
raw <- read.delim(arg("--counts"), check.names=FALSE)
samples <- read.delim(arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
if (!"gene_id" %in% names(raw) || anyDuplicated(raw$gene_id)) stop("Counts need unique gene_id values")
columns <- setdiff(names(raw), "gene_id")
if (anyDuplicated(columns) || anyDuplicated(samples$sample_id)) stop("Duplicate sample IDs")
if (!all(columns %in% samples$sample_id)) stop("Count columns missing from sample metadata")
samples <- samples[match(columns, samples$sample_id), , drop=FALSE]
rownames(samples) <- samples$sample_id
y <- as.matrix(raw[, columns, drop=FALSE]); rownames(y) <- raw$gene_id
if (any(!is.finite(y)) || any(y < 0) || any(y != round(y))) stop("Expected nonnegative integer raw counts")
storage.mode(y) <- "integer"
outdir <- arg("--outdir"); dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
min_count <- as.numeric(arg("--min-count", "5"))
min_total <- as.numeric(arg("--min-total-count", "15"))
fdr <- as.numeric(arg("--fdr", "0.05"))
threshold <- as.numeric(arg("--lfc-threshold", "1"))
write_tsv <- function(frame, target) write.table(frame, target, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
prepare <- function(metadata, text) {
  formula <- as.formula(text)
  for (name in all.vars(formula)) {
    if (!name %in% names(metadata) || any(is.na(metadata[[name]]) | metadata[[name]] == "")) stop(paste("Invalid design column", name))
    metadata[[name]] <- droplevels(factor(metadata[[name]]))
  }
  model <- model.matrix(formula, metadata)
  if (qr(model)$rank < ncol(model) || nrow(model) <= ncol(model)) stop("Design is rank-deficient or has no residual degrees of freedom")
  list(metadata=metadata, formula=formula, model=model)
}
filter_counts <- function(x, design) {
  libs <- colSums(x)
  if (any(!is.finite(libs) | libs <= 0)) stop("A sample has zero/invalid assigned library size")
  keep <- edgeR::filterByExpr(x, design=design, lib.size=libs, min.count=min_count,
                             min.total.count=min_total, large.n=10, min.prop=0.7)
  n <- 1/max(rowSums(qr.Q(qr(design))^2))
  if (n > 10) n <- 10+(n-10)*0.7
  cutoff <- min_count/median(libs)*1e6
  passing <- rowSums(edgeR::cpm(x, lib.size=libs) >= cutoff)
  reason <- ifelse(keep, "eligible", ifelse(rowSums(x) < min_total, "total_count", "insufficient_depth_adjusted_sample_support"))
  list(keep=keep, libs=libs, audit=data.frame(gene_id=rownames(x), total_count=rowSums(x),
       samples_passing=passing, required_samples=ceiling(n-1e-14), cpm_cutoff=cutoff,
       eligible=keep, reason=reason), n=n, cutoff=cutoff)
}

if (arg("--mode", "contrast") == "qc") {
  design <- prepare(samples, arg("--design"))
  filter <- filter_counts(y, design$model)
  if (sum(filter$keep) < 2) stop("Fewer than two genes eligible for expression QC")
  qc <- file.path(outdir, "qc"); dir.create(qc, recursive=TRUE, showWarnings=FALSE)
  all_dds <- DESeqDataSetFromMatrix(y, design$metadata, design$formula)
  all_dds <- estimateSizeFactors(all_dds, type="poscounts")
  norm <- counts(all_dds, normalized=TRUE)
  write_tsv(data.frame(gene_id=rownames(norm), norm, check.names=FALSE), file.path(qc,"normalized_counts.tsv"))
  write_tsv(data.frame(sample_id=columns, assigned_counts=colSums(y), size_factor=sizeFactors(all_dds)), file.path(qc,"library_sizes.tsv"))
  dds <- estimateDispersions(all_dds[filter$keep, ], quiet=TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
  percent <- 100*attr(pca,"percentVar")
  pca$sample_id <- rownames(pca)
  write_tsv(pca, file.path(qc,"pca.tsv"))
  plot <- ggplot(pca,aes(PC1,PC2,color=condition,label=name))+geom_point(size=3)+
    geom_text(vjust=-0.7,check_overlap=TRUE)+theme_bw()+
    labs(title="Extended-exon expression: VST PCA",x=sprintf("PC1 %.1f%%",percent[1]),y=sprintf("PC2 %.1f%%",percent[2]))
  ggsave(file.path(qc,"pca.png"),plot,width=7,height=5,dpi=150)
  ggsave(file.path(qc,"pca.pdf"),plot,width=7,height=5)
  distances <- as.matrix(dist(t(assay(vsd))))
  write_tsv(data.frame(sample_id=rownames(distances),distances,check.names=FALSE),file.path(qc,"sample_distances.tsv"))
  plot <- sample_distance_plot(distances)
  width <- max(8,3+0.28*ncol(distances)); height <- max(7,2+0.28*nrow(distances))
  ggsave(file.path(qc,"sample_distances.png"),plot,width=width,height=height,dpi=150)
  ggsave(file.path(qc,"sample_distances.pdf"),plot,width=width,height=height)
  writeLines(capture.output(sessionInfo()),file.path(qc,"sessionInfo.txt"))
  quit(save="no",status=0)
}

id <- arg("--contrast-id")
con <- contrasts[contrasts$contrast_id == id, , drop=FALSE]
if (nrow(con) != 1) stop("Expected one contrast")
metadata <- samples[samples$condition %in% c(con$numerator,con$denominator), ,drop=FALSE]
design_text <- if ("resolved_design" %in% names(con) && nzchar(con$resolved_design)) con$resolved_design else arg("--design")
design <- prepare(metadata,design_text)
x <- y[,rownames(metadata),drop=FALSE]
filter <- filter_counts(x,design$model)
directory <- file.path(outdir,id); dir.create(directory,recursive=TRUE,showWarnings=FALSE)
write_tsv(filter$audit,file.path(directory,"eligibility.tsv"))
write_tsv(data.frame(sample_id=colnames(x),library_size=filter$libs,
  minimum_raw_count=ceiling(min_count*filter$libs/median(filter$libs)),min_count=min_count,min_total_count=min_total),
  file.path(directory,"filter_parameters.tsv"))
if (sum(filter$keep) < 2) stop(paste("Fewer than two eligible genes in",id))
fit_signature <- list(schema=1,code_hash=arg("--model-code-hash", "standalone"),counts_md5=unname(tools::md5sum(arg("--counts"))),
  sample_metadata=design$metadata,design=design_text,min_count=min_count,min_total=min_total,
  DESeq2=as.character(packageVersion("DESeq2")),edgeR=as.character(packageVersion("edgeR")))
model_file <- file.path(directory,"model.rds")
cached <- if (file.exists(model_file)) tryCatch(readRDS(model_file),error=function(e) NULL) else NULL
if (arg("--force-fit", "false") != "true" && is.list(cached) && identical(cached$signature,fit_signature)) {
  dds <- cached$dds
  message("Reusing one fitted model for both effect hypotheses: ",id)
} else {
  dds <- DESeqDataSetFromMatrix(x[filter$keep,,drop=FALSE],design$metadata,design$formula)
  dds <- DESeq(dds,sfType="poscounts",quiet=TRUE)
  saveRDS(list(signature=fit_signature,dds=dds),model_file)
}
contrast <- c(con$factor,con$numerator,con$denominator)
relaxed <- as.data.frame(results(dds,contrast=contrast,alpha=fdr,independentFiltering=FALSE))
strict <- as.data.frame(results(dds,contrast=contrast,alpha=fdr,lfcThreshold=threshold,
                               altHypothesis="greaterAbs",independentFiltering=FALSE))
valid <- is.finite(relaxed$pvalue) & is.finite(relaxed$stat) & is.finite(strict$pvalue)
for (name in c("relaxed","strict")) {
  table <- get(name)
  table$padj <- NA_real_
  table$padj[valid] <- p.adjust(table$pvalue[valid],method="BH")
  table$gene_id <- rownames(table)
  table$valid_test <- valid
  assign(name,table[,c("gene_id",setdiff(names(table),"gene_id"))])
}
write_tsv(data.frame(gene_id=relaxed$gene_id,valid_test=valid,
  reason=ifelse(valid,"tested","model_or_outlier_test_unavailable")),file.path(directory,"test_universe.tsv"))
ranking <- relaxed
ranking$analysis_mode <- "shared_gsea"
write_tsv(ranking,file.path(directory,"ranking.tsv"))
normalized <- counts(dds,normalized=TRUE)
write_tsv(data.frame(gene_id=rownames(normalized),normalized,check.names=FALSE),file.path(directory,"normalized_counts.tsv"))
writeLines(capture.output(sessionInfo()),file.path(directory,"sessionInfo.txt"))

index <- list()
for (mode in c("stringent","relaxed_fdr","shared_gsea")) {
  if (mode == "stringent" && arg("--stringent","true") != "true") next
  if (mode == "relaxed_fdr" && arg("--relaxed","true") != "true") next
  target <- file.path(directory,mode); dir.create(target,recursive=TRUE,showWarnings=FALSE)
  mode_threshold <- if (mode == "stringent") threshold else 0
  hypothesis <- if (mode == "stringent") sprintf("abs(log2FC) > %g",threshold) else "log2FC != 0"
  table <- if (mode == "stringent") strict else relaxed
  ma_png <- volcano_png <- ma_pdf <- volcano_pdf <- ""
  result_file <- file.path(directory,"ranking.tsv")
  significant <- !is.na(table$padj) & table$padj <= fdr
  if (mode != "shared_gsea") {
    table$analysis_mode <- mode
    table$hypothesis <- hypothesis
    table$significant <- significant
    result_file <- file.path(target,"deseq2.tsv")
    write_tsv(table,result_file)
    write_tsv(table[significant,,drop=FALSE],file.path(target,"significant_genes.tsv"))
    table$direction <- ifelse(significant,ifelse(table$log2FoldChange > 0,"higher_numerator","higher_denominator"),"not_significant")
    table$minus_log10_padj <- -log10(pmax(table$padj,.Machine$double.xmin))
    label <- paste(id,mode,hypothesis,sep=" | ")
    subtitle <- paste("Positive log2FC:",con$numerator,"; negative:",con$denominator,"; BH adjusted p <=",fdr)
    volcano <- ggplot(table,aes(log2FoldChange,minus_log10_padj,color=direction))+geom_point(alpha=.55,size=1,na.rm=TRUE)+
      theme_bw()+labs(title=label,subtitle=subtitle,x="log2 fold change (MLE)",y="-log10 adjusted p")+
      geom_hline(yintercept=-log10(fdr),linetype="dashed")
    ma <- ggplot(table,aes(baseMean,log2FoldChange,color=direction))+geom_point(alpha=.55,size=1,na.rm=TRUE)+
      scale_x_log10()+theme_bw()+labs(title=label,subtitle=subtitle,x="Mean normalized count",y="log2 fold change (MLE)")
    if (mode == "stringent") {
      volcano <- volcano+geom_vline(xintercept=c(-threshold,threshold),linetype="dashed")
      ma <- ma+geom_hline(yintercept=c(-threshold,threshold),linetype="dashed")
    }
    ma_png <- file.path(target,"MA.png"); ma_pdf <- file.path(target,"MA.pdf")
    volcano_png <- file.path(target,"volcano.png"); volcano_pdf <- file.path(target,"volcano.pdf")
    ggsave(ma_png,ma,width=9,height=6,dpi=150); ggsave(ma_pdf,ma,width=9,height=6)
    ggsave(volcano_png,volcano,width=9,height=6,dpi=150); ggsave(volcano_pdf,volcano,width=9,height=6)
  }
  index[[length(index)+1L]] <- data.frame(expression_method="extended_exons_v1",analysis_mode=mode,
    contrast_id=id,numerator=con$numerator,denominator=con$denominator,result_file=result_file,
    tested_genes=sum(valid),eligible_genes=sum(filter$keep),significant=if(mode=="shared_gsea") NA else sum(significant),
    hypothesis=hypothesis,lfc_threshold=mode_threshold,fdr=fdr,
    ma_png=ma_png,ma_pdf=ma_pdf,volcano_png=volcano_png,volcano_pdf=volcano_pdf,
    design_mode=con$design_mode,resolved_design=design_text,paired=con$paired,n_pairs=con$n_pairs,
    warning=con$design_status,stringsAsFactors=FALSE)
}
write_tsv(do.call(rbind,index),arg("--index-file",file.path(directory,"result_index.tsv")))
