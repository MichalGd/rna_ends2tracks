args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name))
  args[[pos + 1]]
}
required <- c("--counts", "--samples", "--contrasts", "--design", "--outdir")
for (item in required) if (is.null(get_arg(item))) stop(paste("Required argument:", item))
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})

count_file <- get_arg("--counts")
samples <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
design_text <- get_arg("--design")
outdir <- get_arg("--outdir")
fdr <- as.numeric(get_arg("--fdr", "0.05"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

raw <- read.delim(count_file, comment.char="#", check.names=FALSE)
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(annotation_columns %in% colnames(raw))) stop("Unexpected featureCounts table")
matrix <- as.matrix(raw[, setdiff(colnames(raw), annotation_columns), drop=FALSE])
rownames(matrix) <- raw$Geneid
storage.mode(matrix) <- "integer"
sample_ids <- samples$sample_id
if (ncol(matrix) != length(sample_ids)) stop("BAM/count column count does not match biological samples")
colnames(matrix) <- sample_ids
rownames(samples) <- samples$sample_id
samples <- samples[colnames(matrix), , drop=FALSE]

design_formula <- as.formula(design_text)
needed <- all.vars(design_formula)
missing <- setdiff(needed, colnames(samples))
if (length(missing)) stop(paste("Design columns missing:", paste(missing, collapse=", ")))
for (name in needed) {
  if (any(is.na(samples[[name]]) | samples[[name]] == "")) stop(paste("Missing design value:", name))
  samples[[name]] <- factor(samples[[name]])
}
model <- model.matrix(design_formula, data=samples)
if (qr(model)$rank < ncol(model)) stop("Design matrix is not full rank; check confounding")
dds <- DESeqDataSetFromMatrix(matrix, samples, design_formula)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)
saveRDS(dds, file.path(outdir, "deseq2_model.rds"))
write.table(counts(dds, normalized=TRUE), file.path(outdir, "normalized_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)

vsd <- vst(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup=intersect(c("condition", "batch"), colnames(samples)), returnData=TRUE)
percent <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) +
  xlab(paste0("PC1: ", percent[1], "%")) + ylab(paste0("PC2: ", percent[2], "%")) + theme_bw()
ggsave(file.path(outdir, "vst_pca.pdf"), p, width=7, height=5)

index <- list()
for (i in seq_len(nrow(contrasts))) {
  contrast <- contrasts[i, ]
  id <- contrast$contrast_id
  result <- results(dds, contrast=c(contrast$factor, contrast$numerator, contrast$denominator), alpha=fdr)
  shrunk <- lfcShrink(dds, contrast=c(contrast$factor, contrast$numerator, contrast$denominator), res=result, type="normal")
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$log2FoldChange_shrunken <- shrunk$log2FoldChange
  table <- table[, c("gene_id", setdiff(colnames(table), "gene_id"))]
  target <- file.path(outdir, paste0(id, ".deseq2.tsv"))
  write.table(table, target, sep="\t", quote=FALSE, row.names=FALSE)
  pdf(file.path(outdir, paste0(id, ".MA.pdf"))); plotMA(result, alpha=fdr); dev.off()
  index[[length(index)+1]] <- data.frame(contrast_id=id, result_file=target,
    significant=sum(!is.na(table$padj) & table$padj < fdr), warning=contrast$design_status)
}
write.table(do.call(rbind, index), file.path(outdir, "result_index.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
