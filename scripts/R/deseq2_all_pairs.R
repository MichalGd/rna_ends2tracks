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
samples_all <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
default_design_text <- get_arg("--design")
outdir <- get_arg("--outdir")
fdr <- as.numeric(get_arg("--fdr", "0.05"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

raw <- read.delim(count_file, comment.char="#", check.names=FALSE)
annotation_columns <- c("Geneid", "Chr", "Start", "End", "Strand", "Length")
if (!all(annotation_columns %in% colnames(raw))) stop("Unexpected featureCounts table")
matrix_all <- as.matrix(raw[, setdiff(colnames(raw), annotation_columns), drop=FALSE])
rownames(matrix_all) <- raw$Geneid
storage.mode(matrix_all) <- "integer"
sample_ids <- samples_all$sample_id
if (ncol(matrix_all) != length(sample_ids)) stop("BAM/count column count does not match biological samples")
colnames(matrix_all) <- sample_ids
rownames(samples_all) <- samples_all$sample_id
samples_all <- samples_all[colnames(matrix_all), , drop=FALSE]

prepare_design <- function(sample_data, design_text, label) {
  design_formula <- as.formula(design_text)
  needed <- all.vars(design_formula)
  missing <- setdiff(needed, colnames(sample_data))
  if (length(missing)) stop(paste("Design columns missing for", label, ":", paste(missing, collapse=", ")))
  for (name in needed) {
    if (any(is.na(sample_data[[name]]) | sample_data[[name]] == "")) {
      stop(paste("Missing design value for", label, ":", name))
    }
    sample_data[[name]] <- factor(sample_data[[name]])
  }
  model <- model.matrix(design_formula, data=sample_data)
  if (qr(model)$rank < ncol(model)) stop(paste("Design matrix is not full rank for", label))
  list(samples=sample_data, formula=design_formula)
}

# The global default-design fit supplies normalized-count and PCA QC only.
global <- prepare_design(samples_all, default_design_text, "global QC")
qc_dds <- DESeqDataSetFromMatrix(matrix_all, global$samples, global$formula)
qc_dds <- qc_dds[rowSums(counts(qc_dds)) >= 10, ]
qc_dds <- DESeq(qc_dds)
saveRDS(qc_dds, file.path(outdir, "deseq2_model.rds"))
write.table(counts(qc_dds, normalized=TRUE), file.path(outdir, "normalized_counts.tsv"),
            sep="\t", quote=FALSE, col.names=NA)

vsd <- vst(qc_dds, blind=FALSE)
pca_groups <- intersect(c("condition", "batch"), colnames(global$samples))
pca <- plotPCA(vsd, intgroup=pca_groups, returnData=TRUE)
percent <- round(100 * attr(pca, "percentVar"))
p <- ggplot(pca, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) +
  xlab(paste0("PC1: ", percent[1], "%")) + ylab(paste0("PC2: ", percent[2], "%")) + theme_bw()
ggsave(file.path(outdir, "vst_pca.pdf"), p, width=7, height=5)

index <- list()
for (i in seq_len(nrow(contrasts))) {
  contrast <- contrasts[i, ]
  id <- contrast$contrast_id
  keep <- samples_all$condition %in% c(contrast$denominator, contrast$numerator)
  sample_data <- samples_all[keep, , drop=FALSE]
  design_text <- default_design_text
  if ("resolved_design" %in% colnames(contrasts) && !is.na(contrast$resolved_design) && nzchar(contrast$resolved_design)) {
    design_text <- contrast$resolved_design
  }
  resolved <- prepare_design(sample_data, design_text, id)
  pair_matrix <- matrix_all[, rownames(resolved$samples), drop=FALSE]
  dds <- DESeqDataSetFromMatrix(pair_matrix, resolved$samples, resolved$formula)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds)
  saveRDS(dds, file.path(outdir, paste0(id, ".deseq2_model.rds")))
  result <- results(dds, contrast=c(contrast$factor, contrast$numerator, contrast$denominator), alpha=fdr)
  shrunk <- lfcShrink(dds, contrast=c(contrast$factor, contrast$numerator, contrast$denominator),
                     res=result, type="normal")
  table <- as.data.frame(result)
  table$gene_id <- rownames(table)
  table$log2FoldChange_shrunken <- shrunk$log2FoldChange
  table <- table[, c("gene_id", setdiff(colnames(table), "gene_id"))]
  target <- file.path(outdir, paste0(id, ".deseq2.tsv"))
  write.table(table, target, sep="\t", quote=FALSE, row.names=FALSE)
  pdf(file.path(outdir, paste0(id, ".MA.pdf"))); plotMA(result, alpha=fdr); dev.off()
  index[[length(index)+1]] <- data.frame(
    contrast_id=id,
    result_file=target,
    significant=sum(!is.na(table$padj) & table$padj < fdr),
    design_mode=if ("design_mode" %in% colnames(contrasts)) contrast$design_mode else "legacy",
    resolved_design=design_text,
    paired=if ("paired" %in% colnames(contrasts)) contrast$paired else FALSE,
    n_pairs=if ("n_pairs" %in% colnames(contrasts)) contrast$n_pairs else 0,
    warning=contrast$design_status,
    check.names=FALSE
  )
}
write.table(do.call(rbind, index), file.path(outdir, "result_index.tsv"),
            sep="\t", quote=FALSE, row.names=FALSE)
