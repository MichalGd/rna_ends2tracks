args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args); if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name)); args[[pos + 1]]
}
suppressPackageStartupMessages({library(DRIMSeq); library(stageR)})
counts_wide <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
catalog <- read.delim(get_arg("--catalog"), check.names=FALSE, stringsAsFactors=FALSE)
samples_all <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
outdir <- get_arg("--outdir"); alpha <- as.numeric(get_arg("--fdr", "0.05"))
default_design_text <- get_arg("--design", "~ condition")
contrast_id <- get_arg("--contrast-id", "")
index_file <- get_arg("--index-file", file.path(outdir, "result_index.tsv"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
rownames(counts_wide) <- counts_wide$pas_id
sample_columns <- intersect(samples_all$sample_id, colnames(counts_wide))
matrix_all <- as.matrix(counts_wide[, sample_columns, drop=FALSE]); storage.mode(matrix_all) <- "integer"
gene_ids <- catalog$gene_id[match(rownames(matrix_all), catalog$pas_id)]
index <- list()
if (nzchar(contrast_id)) {
  contrasts <- contrasts[contrasts$contrast_id == contrast_id, , drop=FALSE]
  if (nrow(contrasts) != 1) stop(paste("Expected exactly one contrast:", contrast_id))
}
for (i in seq_len(nrow(contrasts))) {
  con <- contrasts[i, ]; keep <- samples_all$condition %in% c(con$denominator, con$numerator)
  design_text <- default_design_text
  if ("resolved_design" %in% colnames(contrasts) && !is.na(con$resolved_design) && nzchar(con$resolved_design)) {
    design_text <- con$resolved_design
  }
  design_variables <- all.vars(as.formula(design_text))
  missing <- setdiff(design_variables, colnames(samples_all))
  if (length(missing)) stop(paste("Design columns missing for", con$contrast_id, ":", paste(missing, collapse=", ")))
  samples <- samples_all[keep, unique(c("sample_id", design_variables)), drop=FALSE]
  samples$condition <- relevel(factor(samples$condition), ref=con$denominator)
  for (variable in design_variables) {
    if (any(is.na(samples[[variable]]) | samples[[variable]] == "")) {
      stop(paste("Missing design value for", con$contrast_id, ":", variable))
    }
    if (variable != "condition") samples[[variable]] <- factor(samples[[variable]])
  }
  mat <- matrix_all[, samples$sample_id, drop=FALSE]
  counts <- data.frame(gene_id=gene_ids, feature_id=rownames(mat), mat, check.names=FALSE)
  counts <- counts[counts$gene_id != "" & !is.na(counts$gene_id), , drop=FALSE]
  d <- dmDSdata(counts=counts, samples=samples)
  d <- dmFilter(d, min_samps_feature_expr=2, min_feature_expr=5,
                min_samps_feature_prop=2, min_feature_prop=0.05,
                min_samps_gene_expr=2, min_gene_expr=10)
  design <- model.matrix(as.formula(design_text), data=samples)
  if (qr(design)$rank < ncol(design)) stop(paste("Pair-specific design is not full rank:", con$contrast_id))
  condition_coef <- grep("^condition", colnames(design))
  if (length(condition_coef) != 1) stop(paste("Cannot identify one pairwise condition coefficient:", con$contrast_id))
  d <- dmPrecision(d, design=design); d <- dmFit(d, design=design); d <- dmTest(d, coef=condition_coef)
  gene_result <- DRIMSeq::results(d, level="gene"); feature_result <- DRIMSeq::results(d, level="feature")
  p_screen <- gene_result$pvalue; names(p_screen) <- gene_result$gene_id
  p_confirmation <- matrix(feature_result$pvalue, ncol=1,
                           dimnames=list(feature_result$feature_id, "site"))
  tx2gene <- feature_result[, c("feature_id", "gene_id"), drop=FALSE]
  colnames(tx2gene) <- c("txID", "geneID")
  staged <- stageRTx(pScreen=p_screen, pConfirmation=p_confirmation, pScreenAdjusted=FALSE, tx2gene=tx2gene)
  staged <- stageWiseAdjustment(staged, method="dtu", alpha=alpha)
  adjusted <- as.data.frame(getAdjustedPValues(staged, order=FALSE, onlySignificantGenes=FALSE))
  id_column <- intersect(c("txID", "transcript", "feature_id"), colnames(adjusted))
  adjusted_ids <- if (length(id_column)) adjusted[[id_column[1]]] else rownames(adjusted)
  p_columns <- setdiff(colnames(adjusted), c("txID", "geneID", "transcript", "gene", "feature_id"))
  feature_result$stageR_adjusted <- adjusted[[tail(p_columns, 1)]][match(feature_result$feature_id, adjusted_ids)]
  den <- rowMeans(mat[, samples$condition == con$denominator, drop=FALSE])
  num <- rowMeans(mat[, samples$condition == con$numerator, drop=FALSE])
  den_total <- ave(den, gene_ids, FUN=sum); num_total <- ave(num, gene_ids, FUN=sum)
  feature_result$PAU_denominator <- den[match(feature_result$feature_id, rownames(mat))] / pmax(den_total[match(feature_result$feature_id, rownames(mat))], 1e-12)
  feature_result$PAU_numerator <- num[match(feature_result$feature_id, rownames(mat))] / pmax(num_total[match(feature_result$feature_id, rownames(mat))], 1e-12)
  feature_result$delta_PAU <- feature_result$PAU_numerator - feature_result$PAU_denominator
  positions <- catalog$start[match(feature_result$feature_id, catalog$pas_id)]
  result_genes <- feature_result$gene_id
  feature_result$weighted_genomic_position_shift_nt <- ave(feature_result$PAU_numerator * positions, result_genes, FUN=sum) - ave(feature_result$PAU_denominator * positions, result_genes, FUN=sum)
  target <- file.path(outdir, paste0(con$contrast_id, ".drimseq_stager.tsv"))
  write.table(feature_result, target, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(gene_result, file.path(outdir, paste0(con$contrast_id, ".gene_screen.tsv")), sep="\t", quote=FALSE, row.names=FALSE)
  index[[length(index)+1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=target, tested_sites=nrow(feature_result),
    confirmed_sites=sum(!is.na(feature_result$stageR_adjusted) & feature_result$stageR_adjusted < alpha),
    design_mode=if ("design_mode" %in% colnames(contrasts)) con$design_mode else "legacy",
    resolved_design=design_text,
    paired=if ("paired" %in% colnames(contrasts)) con$paired else FALSE,
    n_pairs=if ("n_pairs" %in% colnames(contrasts)) con$n_pairs else 0,
    warning=con$design_status, check.names=FALSE
  )
}
write.table(do.call(rbind, index), index_file, sep="\t", quote=FALSE, row.names=FALSE)
