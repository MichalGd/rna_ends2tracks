args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args); if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name)); args[[pos + 1]]
}
suppressPackageStartupMessages(library(DEXSeq))
counts_df <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
catalog <- read.delim(get_arg("--catalog"), check.names=FALSE, stringsAsFactors=FALSE)
samples <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
outdir <- get_arg("--outdir"); min_count <- as.integer(get_arg("--min-count", "5"))
default_design_text <- get_arg("--design", "~ condition")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
meta <- merge(catalog[, c("pas_id", "gene_id")], counts_df[, "pas_id", drop=FALSE], by="pas_id", sort=FALSE)
rownames(counts_df) <- counts_df$pas_id
count_matrix <- as.matrix(counts_df[, setdiff(colnames(counts_df), "pas_id"), drop=FALSE]); storage.mode(count_matrix) <- "integer"
index <- list()
for (i in seq_len(nrow(contrasts))) {
  con <- contrasts[i, ]; keep_samples <- samples$sample_id[samples$condition %in% c(con$denominator, con$numerator)]
  sample_data <- samples[match(keep_samples, samples$sample_id), , drop=FALSE]
  rownames(sample_data) <- sample_data$sample_id
  sample_data$condition <- relevel(factor(sample_data$condition), ref=con$denominator)
  design_text <- default_design_text
  if ("resolved_design" %in% colnames(contrasts) && !is.na(con$resolved_design) && nzchar(con$resolved_design)) {
    design_text <- con$resolved_design
  }
  design_variables <- all.vars(as.formula(design_text))
  covariates <- setdiff(design_variables, "condition")
  missing <- setdiff(design_variables, colnames(sample_data))
  if (length(missing)) stop(paste("Design columns missing for", con$contrast_id, ":", paste(missing, collapse=", ")))
  for (variable in design_variables) {
    if (any(is.na(sample_data[[variable]]) | sample_data[[variable]] == "")) {
      stop(paste("Missing design value for", con$contrast_id, ":", variable))
    }
    if (variable != "condition") sample_data[[variable]] <- factor(sample_data[[variable]])
  }
  pair_design <- model.matrix(as.formula(design_text), data=sample_data)
  if (qr(pair_design)$rank < ncol(pair_design)) stop(paste("Pair-specific design is not full rank:", con$contrast_id))
  covariate_terms <- if (length(covariates)) paste0(covariates, ":exon") else character(0)
  full_formula <- as.formula(paste("~ sample + exon", paste(c(covariate_terms, "condition:exon"), collapse=" + "), sep=" + "))
  reduced_formula <- as.formula(paste("~ sample + exon", paste(covariate_terms, collapse=" + "), sep=if (length(covariate_terms)) " + " else ""))
  keep_sites <- catalog$confidence %in% c("high_confidence", "rescued_a_rich") & catalog$gene_id != ""
  gene_sites <- table(catalog$gene_id[keep_sites]); eligible_genes <- names(gene_sites[gene_sites >= 2])
  keep_sites <- keep_sites & catalog$gene_id %in% eligible_genes & rowSums(count_matrix[, keep_samples, drop=FALSE]) >= min_count
  gene_sites <- table(catalog$gene_id[keep_sites]); eligible_genes <- names(gene_sites[gene_sites >= 2])
  keep_sites <- keep_sites & catalog$gene_id %in% eligible_genes
  selected <- catalog$pas_id[keep_sites]
  if (length(selected) < 2) stop(paste("No testable APA-A sites for", con$contrast_id))
  mat <- count_matrix[selected, keep_samples, drop=FALSE]
  dxd <- DEXSeqDataSet(mat, sampleData=sample_data, design=full_formula,
                       featureID=selected, groupID=catalog$gene_id[match(selected, catalog$pas_id)])
  dxd <- estimateSizeFactors(dxd); dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd, reducedModel=reduced_formula)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
  result <- as.data.frame(DEXSeqResults(dxd)); result$pas_id <- result$featureID
  norm <- counts(dxd, normalized=TRUE)
  den <- rowMeans(norm[, sample_data$condition == con$denominator, drop=FALSE])
  num <- rowMeans(norm[, sample_data$condition == con$numerator, drop=FALSE])
  groups <- catalog$gene_id[match(rownames(norm), catalog$pas_id)]
  den_total <- ave(den, groups, FUN=sum); num_total <- ave(num, groups, FUN=sum)
  result$PAU_denominator <- den / pmax(den_total, 1e-12)
  result$PAU_numerator <- num / pmax(num_total, 1e-12)
  result$delta_PAU <- result$PAU_numerator - result$PAU_denominator
  positions <- catalog$start[match(result$pas_id, catalog$pas_id)]
  result$weighted_genomic_position_shift_nt <- ave(result$PAU_numerator * positions, groups, FUN=sum) - ave(result$PAU_denominator * positions, groups, FUN=sum)
  result <- result[, c("pas_id", setdiff(colnames(result), "pas_id"))]
  target <- file.path(outdir, paste0(con$contrast_id, ".dexseq.tsv"))
  write.table(result, target, sep="\t", quote=FALSE, row.names=FALSE)
  index[[length(index)+1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=target, tested_sites=nrow(result),
    significant_sites=sum(!is.na(result$padj) & result$padj < 0.05),
    design_mode=if ("design_mode" %in% colnames(contrasts)) con$design_mode else "legacy",
    resolved_design=design_text,
    paired=if ("paired" %in% colnames(contrasts)) con$paired else FALSE,
    n_pairs=if ("n_pairs" %in% colnames(contrasts)) con$n_pairs else 0,
    warning=con$design_status, check.names=FALSE
  )
}
write.table(do.call(rbind, index), file.path(outdir, "result_index.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
