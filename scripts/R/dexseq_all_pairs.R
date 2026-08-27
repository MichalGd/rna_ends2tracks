args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name))
  args[[pos + 1]]
}

dexseq_formulas <- function(design_text) {
  variables <- all.vars(as.formula(design_text))
  covariates <- setdiff(variables, "condition")
  covariate_interactions <- if (length(covariates)) paste0(covariates, ":exon") else character(0)
  full_terms <- c(covariate_interactions, "condition:exon")
  full_formula <- as.formula(paste("~ sample + exon +", paste(full_terms, collapse=" + ")))
  reduced_formula <- if (length(covariate_interactions)) {
    as.formula(paste("~ sample + exon +", paste(covariate_interactions, collapse=" + ")))
  } else as.formula("~ sample + exon")
  list(
    variables=variables,
    covariates=covariates,
    full=full_formula,
    reduced=reduced_formula
  )
}

tabular_dexseq_result <- function(result) {
  list_columns <- vapply(result, is.list, logical(1))
  if (any(list_columns)) {
    message(
      "Dropping non-tabular DEXSeq result columns: ",
      paste(colnames(result)[list_columns], collapse=", ")
    )
    result <- result[, !list_columns, drop=FALSE]
  }
  result
}

pooled_raw_for_pas <- function(pooled_counts, pas_ids) {
  ids <- as.character(pas_ids)
  missing <- setdiff(ids, names(pooled_counts))
  if (length(missing)) {
    stop(
      "Comparator PAS absent from the C3 count matrix: ",
      paste(head(missing, 10), collapse=", ")
    )
  }
  unname(pooled_counts[match(ids, names(pooled_counts))])
}

catalog_rows_for_pas <- function(catalog, pas_ids) {
  ids <- as.character(pas_ids)
  positions <- match(ids, as.character(catalog$pas_id))
  if (anyNA(positions)) {
    missing <- unique(ids[is.na(positions)])
    stop(
      "DEXSeq PAS absent from the active-PAS catalog: ",
      paste(head(missing, 10), collapse=", ")
    )
  }
  positions
}

name_count_rows_by_pas <- function(count_matrix, pas_ids) {
  ids <- as.character(pas_ids)
  if (nrow(count_matrix) != length(ids)) {
    stop("Normalized DEXSeq count rows do not match the selected PAS universe")
  }
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    stop("Selected PAS identifiers must be complete and unique")
  }
  rownames(count_matrix) <- ids
  count_matrix
}

if ("--self-test" %in% args) {
  unpaired <- dexseq_formulas("~ condition")
  stopifnot(
    paste(deparse(unpaired$full), collapse=" ") == "~sample + exon + condition:exon",
    paste(deparse(unpaired$reduced), collapse=" ") == "~sample + exon"
  )
  paired <- dexseq_formulas("~ subject + condition")
  stopifnot(
    paste(deparse(paired$full), collapse=" ") == "~sample + exon + subject:exon + condition:exon",
    paste(deparse(paired$reduced), collapse=" ") == "~sample + exon + subject:exon"
  )
  fixture <- data.frame(featureID=c("p1", "p2"), padj=c(0.01, 0.20), stringsAsFactors=FALSE)
  fixture$genomicData <- I(list(list(chr="chr1"), list(chr="chr2")))
  clean <- tabular_dexseq_result(fixture)
  stopifnot(identical(colnames(clean), c("featureID", "padj")))
  target <- tempfile(fileext=".tsv")
  write.table(clean, target, sep="\t", quote=FALSE, row.names=FALSE)
  stopifnot(file.exists(target), nrow(read.delim(target, check.names=FALSE)) == 2)
  unlink(target)
  raw_fixture <- matrix(
    c(1L, 2L, 10L, 20L), nrow=2,
    dimnames=list(c("p1", "p2"), c("s1", "s2"))
  )
  pooled_fixture <- rowSums(raw_fixture)
  factor_ids <- factor(c("p2", "p1"), levels=c("p1", "p2", "unused"))
  stopifnot(identical(pooled_raw_for_pas(pooled_fixture, factor_ids), c(22, 11)))
  catalog_fixture <- data.frame(
    pas_id=c("p1", "p2"), gene_id=c("g1", "g2"), stringsAsFactors=FALSE
  )
  stopifnot(identical(catalog_rows_for_pas(catalog_fixture, factor_ids), c(2L, 1L)))
  dexseq_named_fixture <- raw_fixture
  rownames(dexseq_named_fixture) <- c("g1:p1", "g2:p2")
  restored_fixture <- name_count_rows_by_pas(dexseq_named_fixture, c("p1", "p2"))
  stopifnot(identical(rownames(restored_fixture), c("p1", "p2")))
  missing_error <- tryCatch(
    { pooled_raw_for_pas(pooled_fixture, "absent"); "" },
    error=function(condition) conditionMessage(condition)
  )
  stopifnot(grepl("Comparator PAS absent from the C3 count matrix: absent", missing_error, fixed=TRUE))
  cat("DEXSeq formula, serialization, PAS naming, and comparator lookup self-test: PASS\n")
  quit(save="no", status=0)
}

suppressPackageStartupMessages(library(DEXSeq))

counts_df <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
catalog <- read.delim(get_arg("--catalog"), check.names=FALSE, stringsAsFactors=FALSE)
samples_all <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
contrasts <- read.delim(get_arg("--contrasts"), check.names=FALSE, stringsAsFactors=FALSE)
outdir <- get_arg("--outdir")
min_count <- as.integer(get_arg("--min-count", "5"))
fdr <- as.numeric(get_arg("--fdr", "0.05"))
default_design_text <- get_arg("--design", "~ condition")
contrast_id <- get_arg("--contrast-id", "")
index_file <- get_arg("--index-file", file.path(outdir, "result_index.tsv"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

required_catalog <- c("pas_id", "gene_id", "strand", "summit_start", "assignment_status")
if (!all(required_catalog %in% colnames(catalog))) stop("Active-PAS catalog contract is incomplete")
rownames(counts_df) <- counts_df$pas_id
sample_columns <- intersect(samples_all$sample_id, setdiff(colnames(counts_df), "pas_id"))
count_matrix <- as.matrix(counts_df[, sample_columns, drop=FALSE])
storage.mode(count_matrix) <- "integer"
if (nzchar(contrast_id)) {
  contrasts <- contrasts[contrasts$contrast_id == contrast_id, , drop=FALSE]
  if (nrow(contrasts) != 1) stop(paste("Expected exactly one contrast:", contrast_id))
}

ratio_text <- function(numerator, denominator) {
  if (numerator == 0 && denominator == 0) return("NA")
  if (denominator == 0) return("Inf")
  format(numerator / denominator, digits=12)
}
direction_from_cross_product <- function(dt, pt, dc, pc) {
  if ((dt == 0 && pt == 0) || (dc == 0 && pc == 0)) return("not_classifiable")
  left <- dt * pc; right <- pt * dc
  if (left > right) "distal" else if (left < right) "proximal" else "no_directional_change"
}

index <- list()
for (i in seq_len(nrow(contrasts))) {
  con <- contrasts[i, ]
  keep <- samples_all$sample_id %in% sample_columns & samples_all$condition %in% c(con$denominator, con$numerator)
  sample_data <- samples_all[keep, , drop=FALSE]
  rownames(sample_data) <- sample_data$sample_id
  sample_data$condition <- relevel(factor(sample_data$condition), ref=con$denominator)
  design_text <- if ("resolved_design" %in% colnames(contrasts) && nzchar(con$resolved_design)) con$resolved_design else default_design_text
  formulas <- dexseq_formulas(design_text)
  variables <- formulas$variables
  missing <- setdiff(variables, colnames(sample_data))
  if (length(missing)) stop(paste("Design columns missing for", con$contrast_id, paste(missing, collapse=", ")))
  for (variable in variables) sample_data[[variable]] <- factor(sample_data[[variable]])
  design_matrix <- model.matrix(as.formula(design_text), data=sample_data)
  if (qr(design_matrix)$rank < ncol(design_matrix)) stop(paste("Pair-specific design is not full rank:", con$contrast_id))
  full_formula <- formulas$full
  reduced_formula <- formulas$reduced

  eligible <- catalog$assignment_status == "unique" & !is.na(catalog$gene_id) & nzchar(catalog$gene_id)
  first_count <- table(catalog$gene_id[eligible])
  eligible <- eligible & catalog$gene_id %in% names(first_count[first_count >= 2])
  eligible <- eligible & rowSums(count_matrix[catalog$pas_id, rownames(sample_data), drop=FALSE]) >= min_count
  second_count <- table(catalog$gene_id[eligible])
  eligible <- eligible & catalog$gene_id %in% names(second_count[second_count >= 2])
  selected <- catalog$pas_id[eligible]
  if (length(selected) < 2) stop(paste("No testable APA-A sites for", con$contrast_id))
  mat <- count_matrix[selected, rownames(sample_data), drop=FALSE]
  group_ids <- catalog$gene_id[match(selected, catalog$pas_id)]
  dxd <- DEXSeqDataSet(mat, sampleData=sample_data, design=full_formula, featureID=selected, groupID=group_ids)
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd, reducedModel=reduced_formula)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
  result <- tabular_dexseq_result(as.data.frame(DEXSeqResults(dxd)))
  result$pas_id <- as.character(result$featureID)
  result_catalog_rows <- catalog_rows_for_pas(catalog, result$pas_id)
  norm <- counts(dxd, normalized=TRUE)
  norm <- name_count_rows_by_pas(norm, selected)
  den <- rowMeans(norm[, sample_data$condition == con$denominator, drop=FALSE])
  num <- rowMeans(norm[, sample_data$condition == con$numerator, drop=FALSE])
  all_mean <- rowMeans(norm)
  pooled_raw_counts <- rowSums(count_matrix[, rownames(sample_data), drop=FALSE])
  genes <- catalog$gene_id[match(rownames(norm), catalog$pas_id)]
  den_total <- ave(den, genes, FUN=sum); num_total <- ave(num, genes, FUN=sum)
  result$mean_normalized_count <- all_mean[match(result$pas_id, rownames(norm))]
  result$PAU_denominator <- ifelse(den_total > 0, den / den_total, 0)
  result$PAU_numerator <- ifelse(num_total > 0, num / num_total, 0)
  result$delta_PAU <- result$PAU_numerator - result$PAU_denominator
  result$gene_id <- as.character(catalog$gene_id[result_catalog_rows])
  result$summit_start <- catalog$summit_start[result_catalog_rows]
  result$strand <- as.character(catalog$strand[result_catalog_rows])
  if (anyNA(result$gene_id) || any(!nzchar(result$gene_id))) {
    stop("Active-PAS catalog returned an empty gene_id for a DEXSeq PAS")
  }
  result <- result[, c("pas_id", "gene_id", "summit_start", "strand", setdiff(colnames(result), c("pas_id", "gene_id", "summit_start", "strand")))]
  target <- file.path(outdir, paste0(con$contrast_id, ".dexseq.tsv"))
  write.table(result, target, sep="\t", quote=FALSE, row.names=FALSE)

  shift_rows <- list()
  for (gene in unique(result$gene_id)) {
    rows <- result[result$gene_id == gene, , drop=FALSE]
    significant <- rows[!is.na(rows$padj) & rows$padj < fdr, , drop=FALSE]
    if (!nrow(significant)) {
      shift_rows[[length(shift_rows) + 1]] <- data.frame(gene_id=gene, shift="no_shift", proximal_pas="", distal_pas="",
        ratio_numerator="NA", ratio_denominator="NA", comparator_rule="no_significant_pas", stringsAsFactors=FALSE)
      next
    }
    ordered <- significant[order(significant$padj, -abs(significant$delta_PAU), -significant$mean_normalized_count, significant$summit_start), , drop=FALSE]
    if (nrow(ordered) >= 2) {
      chosen <- ordered[1:2, , drop=FALSE]; comparator_rule <- "two_significant_lowest_padj"
    } else {
      eligible_other <- (
        !is.na(rows$pas_id) & nzchar(as.character(rows$pas_id)) &
        as.character(rows$pas_id) != as.character(ordered$pas_id[1]) &
        !is.na(rows$mean_normalized_count) & is.finite(rows$mean_normalized_count) &
        rows$mean_normalized_count > 0
      )
      other <- rows[eligible_other, , drop=FALSE]
      if (!nrow(other)) {
        shift_rows[[length(shift_rows) + 1]] <- data.frame(gene_id=gene, shift="not_classifiable", proximal_pas="", distal_pas="",
          ratio_numerator="NA", ratio_denominator="NA", comparator_rule="no_nonzero_comparator", stringsAsFactors=FALSE)
        next
      }
      pooled_raw <- pooled_raw_for_pas(pooled_raw_counts, other$pas_id)
      other <- other[order(-other$mean_normalized_count, -pooled_raw, other$summit_start), , drop=FALSE]
      chosen <- rbind(ordered[1, , drop=FALSE], other[1, , drop=FALSE])
      comparator_rule <- "one_significant_vs_highest_mean_normalized_other"
    }
    if (chosen$strand[1] == "+") chosen <- chosen[order(chosen$summit_start), , drop=FALSE]
    else chosen <- chosen[order(-chosen$summit_start), , drop=FALSE]
    proximal <- chosen[1, ]; distal <- chosen[2, ]
    shift_rows[[length(shift_rows) + 1]] <- data.frame(
      gene_id=gene,
      shift=direction_from_cross_product(distal$PAU_numerator, proximal$PAU_numerator, distal$PAU_denominator, proximal$PAU_denominator),
      proximal_pas=proximal$pas_id, distal_pas=distal$pas_id,
      ratio_numerator=ratio_text(distal$PAU_numerator, proximal$PAU_numerator),
      ratio_denominator=ratio_text(distal$PAU_denominator, proximal$PAU_denominator),
      comparator_rule=comparator_rule, stringsAsFactors=FALSE
    )
  }
  shift_target <- file.path(outdir, paste0(con$contrast_id, ".apa_shift.tsv"))
  write.table(do.call(rbind, shift_rows), shift_target, sep="\t", quote=FALSE, row.names=FALSE)
  index[[length(index) + 1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=target, shift_file=shift_target, tested_sites=nrow(result),
    significant_sites=sum(!is.na(result$padj) & result$padj < fdr), design_mode=con$design_mode,
    resolved_design=design_text, paired=con$paired, n_pairs=con$n_pairs, warning=con$design_status,
    check.names=FALSE
  )
}
write.table(do.call(rbind, index), index_file, sep="\t", quote=FALSE, row.names=FALSE)
