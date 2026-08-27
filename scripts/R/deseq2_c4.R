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

raw <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
samples_all <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
default_design_text <- get_arg("--design")
outdir <- get_arg("--outdir")
fdr <- as.numeric(get_arg("--fdr", "0.05"))
plot_lfc <- as.numeric(get_arg("--plot-lfc", "1"))
mode <- get_arg("--mode", "all")
contrast_id <- get_arg("--contrast-id", "")
index_file <- get_arg("--index-file", file.path(outdir, "result_index.tsv"))
factor_output <- get_arg("--factor-output", file.path(outdir, "track_size_factors.tsv"))
if (!mode %in% c("all", "qc", "contrast")) stop("--mode must be all, qc, or contrast")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
if (!"gene_id" %in% colnames(raw)) stop("C4 matrix requires gene_id")
sample_columns <- intersect(samples_all$sample_id, setdiff(colnames(raw), "gene_id"))
matrix_all <- as.matrix(raw[, sample_columns, drop=FALSE])
rownames(matrix_all) <- raw$gene_id
storage.mode(matrix_all) <- "integer"
samples_all <- samples_all[match(sample_columns, samples_all$sample_id), , drop=FALSE]
rownames(samples_all) <- samples_all$sample_id

prepare_design <- function(sample_data, design_text, label) {
  formula <- as.formula(design_text)
  needed <- all.vars(formula)
  missing <- setdiff(needed, colnames(sample_data))
  if (length(missing)) stop(paste("Design columns missing for", label, paste(missing, collapse=", ")))
  for (name in needed) {
    if (any(is.na(sample_data[[name]]) | sample_data[[name]] == "")) stop(paste("Missing", name, "for", label))
    sample_data[[name]] <- factor(sample_data[[name]])
  }
  model <- model.matrix(formula, data=sample_data)
  if (qr(model)$rank < ncol(model)) stop(paste("Design matrix is not full rank for", label))
  list(samples=sample_data, formula=formula)
}

if (mode %in% c("all", "qc")) {
  global <- prepare_design(samples_all, default_design_text, "global C4 QC")
  dds_all <- DESeqDataSetFromMatrix(matrix_all, global$samples, global$formula)
  dds_all <- estimateSizeFactors(dds_all, type="poscounts")
  keep <- rowSums(counts(dds_all)) >= 10
  if (sum(keep) < 2) stop("C4 matrix has fewer than two genes with at least 10 counts")
  dds <- dds_all[keep, ]
  sizeFactors(dds) <- sizeFactors(dds_all)
  dds <- DESeq(dds)
  saveRDS(dds, file.path(outdir, "C4_deseq2_model.rds"))
  normalized <- counts(dds, normalized=TRUE)
  write.table(normalized, file.path(outdir, "C4_normalized_counts.tsv"), sep="\t", quote=FALSE, col.names=NA)
  size_factor <- sizeFactors(dds_all)
  column_sum <- colSums(matrix_all[, names(size_factor), drop=FALSE])
  positive <- column_sum[column_sum > 0]
  geometric_mean <- if (length(positive)) exp(mean(log(positive))) else NA_real_
  factors <- data.frame(
    sample_id=names(size_factor), size_factor=as.numeric(size_factor), C4_column_sum=as.numeric(column_sum),
    cohort_geometric_mean_column_sum=geometric_mean,
    deseq2_scale=1 / as.numeric(size_factor),
    robust_effective_library=as.numeric(size_factor) * geometric_mean,
    robust_cpm_scale=1000000 / (as.numeric(size_factor) * geometric_mean),
    estimator="DESeq2_poscounts", stringsAsFactors=FALSE
  )
  dir.create(dirname(factor_output), recursive=TRUE, showWarnings=FALSE)
  write.table(factors, factor_output, sep="\t", quote=FALSE, row.names=FALSE)
  robust_expected <- DESeq2::fpm(dds_all, robust=TRUE)
  robust_scaled <- sweep(counts(dds_all), 2, factors$robust_effective_library, "/") * 1000000
  tolerance <- 1e-8 * max(1, max(abs(robust_expected)))
  if (max(abs(robust_expected - robust_scaled)) > tolerance) stop("Robust CPM formula does not match DESeq2::fpm(..., robust=TRUE)")
  vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
  pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
  percent <- round(100 * attr(pca, "percentVar"))
  pca$sample_id <- rownames(pca)
  write.table(pca, file.path(outdir, "C4_vst_pca.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
  plot <- ggplot(pca, aes(PC1, PC2, color=condition, label=name)) + geom_point(size=3) +
    geom_text(vjust=-0.7, check_overlap=TRUE) +
    xlab(paste0("PC1: ", percent[1], "%")) + ylab(paste0("PC2: ", percent[2], "%")) + theme_bw()
  ggsave(file.path(outdir, "C4_vst_pca.pdf"), plot, width=7, height=5)
  ggsave(file.path(outdir, "C4_vst_pca.png"), plot, width=7, height=5, dpi=150)
  distance_matrix <- as.matrix(dist(t(assay(vsd))))
  write.table(distance_matrix, file.path(outdir, "C4_sample_distances.tsv"), sep="\t", quote=FALSE, col.names=NA)
  distance_long <- data.frame(
    sample_x=rep(colnames(distance_matrix), each=nrow(distance_matrix)),
    sample_y=rep(rownames(distance_matrix), times=ncol(distance_matrix)),
    distance=as.vector(distance_matrix), stringsAsFactors=FALSE
  )
  distance_plot <- ggplot(distance_long, aes(sample_x, sample_y, fill=distance)) +
    geom_tile() + scale_fill_gradient(low="white", high="#2166AC") +
    labs(x=NULL, y=NULL, fill="VST distance") + theme_bw() +
    theme(axis.text.x=element_text(angle=60, hjust=1))
  ggsave(file.path(outdir, "C4_sample_distances.pdf"), distance_plot, width=8, height=7)
  ggsave(file.path(outdir, "C4_sample_distances.png"), distance_plot, width=8, height=7, dpi=150)
}
if (mode == "qc") quit(save="no", status=0)
if (nzchar(contrast_id)) {
  contrasts <- contrasts[contrasts$contrast_id == contrast_id, , drop=FALSE]
  if (nrow(contrasts) != 1) stop(paste("Expected exactly one contrast:", contrast_id))
}

index <- list()
for (i in seq_len(nrow(contrasts))) {
  con <- contrasts[i, ]
  keep <- samples_all$condition %in% c(con$denominator, con$numerator)
  sample_data <- samples_all[keep, , drop=FALSE]
  design_text <- if ("resolved_design" %in% colnames(contrasts) && nzchar(con$resolved_design)) con$resolved_design else default_design_text
  resolved <- prepare_design(sample_data, design_text, con$contrast_id)
  pair_matrix <- matrix_all[, rownames(resolved$samples), drop=FALSE]
  dds <- DESeqDataSetFromMatrix(pair_matrix, resolved$samples, resolved$formula)
  dds <- dds[rowSums(counts(dds)) >= 10, ]
  dds <- DESeq(dds, sfType="poscounts")
  result <- results(dds, contrast=c(con$factor, con$numerator, con$denominator), alpha=fdr)
  table <- as.data.frame(result); table$gene_id <- rownames(table)
  table <- table[, c("gene_id", setdiff(colnames(table), "gene_id"))]
  target <- file.path(outdir, paste0(con$contrast_id, ".deseq2.tsv"))
  model_target <- file.path(outdir, paste0(con$contrast_id, ".deseq2_model.rds"))
  ma_target <- file.path(outdir, paste0(con$contrast_id, ".MA.pdf"))
  ma_png_target <- file.path(outdir, paste0(con$contrast_id, ".MA.png"))
  volcano_target <- file.path(outdir, paste0(con$contrast_id, ".volcano.pdf"))
  volcano_png_target <- file.path(outdir, paste0(con$contrast_id, ".volcano.png"))
  write.table(table, target, sep="\t", quote=FALSE, row.names=FALSE)
  saveRDS(dds, model_target)
  pdf(ma_target); plotMA(result, alpha=fdr); dev.off()
  png(ma_png_target, width=1050, height=750, res=150); plotMA(result, alpha=fdr); dev.off()
  volcano <- table
  volcano$minus_log10_padj <- -log10(pmax(volcano$padj, .Machine$double.xmin))
  volcano$significance <- ifelse(
    !is.na(volcano$padj) & volcano$padj < fdr & volcano$log2FoldChange >= plot_lfc, "up",
    ifelse(!is.na(volcano$padj) & volcano$padj < fdr & volcano$log2FoldChange <= -plot_lfc, "down", "not_significant")
  )
  volcano_plot <- ggplot(volcano, aes(log2FoldChange, minus_log10_padj, color=significance)) +
    geom_point(alpha=0.55, size=1) +
    scale_color_manual(values=c(down="#2166AC", not_significant="#BDBDBD", up="#B2182B")) +
    geom_vline(xintercept=c(-plot_lfc, plot_lfc), linetype="dashed", color="#555555") +
    geom_hline(yintercept=-log10(fdr), linetype="dashed", color="#555555") +
    labs(title=con$contrast_id, x="log2 fold change", y="-log10 adjusted P", color=NULL) + theme_bw()
  ggsave(volcano_target, volcano_plot, width=7, height=5)
  ggsave(volcano_png_target, volcano_plot, width=7, height=5, dpi=150)
  index[[length(index) + 1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=target, significant=sum(!is.na(table$padj) & table$padj < fdr),
    ma_pdf=ma_target, ma_png=ma_png_target, volcano_pdf=volcano_target, volcano_png=volcano_png_target,
    design_mode=con$design_mode, resolved_design=design_text, paired=con$paired, n_pairs=con$n_pairs,
    warning=con$design_status, stringsAsFactors=FALSE
  )
}
write.table(do.call(rbind, index), index_file, sep="\t", quote=FALSE, row.names=FALSE)
