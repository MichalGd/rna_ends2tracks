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
  list(
    variables=variables,
    full=as.formula(paste("~ sample + exon +", paste(full_terms, collapse=" + "))),
    reduced=if (length(covariate_interactions)) {
      as.formula(paste("~ sample + exon +", paste(covariate_interactions, collapse=" + ")))
    } else as.formula("~ sample + exon")
  )
}

tabular_dexseq_result <- function(result) {
  list_columns <- vapply(result, is.list, logical(1))
  if (any(list_columns)) {
    message("Dropping non-tabular DEXSeq result columns: ",
            paste(colnames(result)[list_columns], collapse=", "))
    result <- result[, !list_columns, drop=FALSE]
  }
  result
}

catalog_rows_for_pas <- function(catalog, pas_ids) {
  positions <- match(as.character(pas_ids), as.character(catalog$pas_id))
  if (anyNA(positions)) {
    missing <- unique(as.character(pas_ids)[is.na(positions)])
    stop("DEXSeq PAS absent from the active-PAS catalog: ", paste(head(missing, 10), collapse=", "))
  }
  positions
}

within_gene_pau <- function(count_matrix, gene_ids) {
  pau <- matrix(
    NA_real_, nrow=nrow(count_matrix), ncol=ncol(count_matrix),
    dimnames=dimnames(count_matrix)
  )
  for (gene in unique(gene_ids)) {
    rows <- which(gene_ids == gene)
    totals <- colSums(count_matrix[rows, , drop=FALSE])
    valid <- is.finite(totals) & totals > 0
    if (any(valid)) {
      pau[rows, valid] <- sweep(count_matrix[rows, valid, drop=FALSE], 2, totals[valid], "/")
    }
  }
  pau
}

row_mean_or_na <- function(values) {
  result <- rowMeans(values, na.rm=TRUE)
  result[rowSums(!is.na(values)) == 0] <- NA_real_
  result
}

unpaired_effect <- function(pau, sample_data, denominator, numerator) {
  den <- pau[, sample_data$condition == denominator, drop=FALSE]
  num <- pau[, sample_data$condition == numerator, drop=FALSE]
  list(
    denominator=row_mean_or_na(den),
    numerator=row_mean_or_na(num),
    delta=row_mean_or_na(num) - row_mean_or_na(den),
    pair_deltas=matrix(numeric(0), nrow=nrow(pau), ncol=0),
    pair_ids=character(0)
  )
}

paired_effect <- function(pau, sample_data, denominator, numerator, pairing_column) {
  den_data <- sample_data[sample_data$condition == denominator, , drop=FALSE]
  num_data <- sample_data[sample_data$condition == numerator, , drop=FALSE]
  pair_ids <- intersect(as.character(den_data[[pairing_column]]), as.character(num_data[[pairing_column]]))
  pair_ids <- sort(unique(pair_ids[nzchar(pair_ids)]))
  if (!length(pair_ids)) stop("Paired APA-A2 contrast has no complete pairs")
  den_columns <- vapply(pair_ids, function(pair_id) {
    hits <- rownames(den_data)[as.character(den_data[[pairing_column]]) == pair_id]
    if (length(hits) != 1) stop("APA-A2 requires one denominator sample per pair: ", pair_id)
    hits[[1]]
  }, character(1))
  num_columns <- vapply(pair_ids, function(pair_id) {
    hits <- rownames(num_data)[as.character(num_data[[pairing_column]]) == pair_id]
    if (length(hits) != 1) stop("APA-A2 requires one numerator sample per pair: ", pair_id)
    hits[[1]]
  }, character(1))
  den <- pau[, den_columns, drop=FALSE]
  num <- pau[, num_columns, drop=FALSE]
  valid <- !is.na(den) & !is.na(num)
  pair_deltas <- num - den
  pair_deltas[!valid] <- NA_real_
  den[!valid] <- NA_real_
  num[!valid] <- NA_real_
  colnames(pair_deltas) <- pair_ids
  list(
    denominator=row_mean_or_na(den),
    numerator=row_mean_or_na(num),
    delta=row_mean_or_na(pair_deltas),
    pair_deltas=pair_deltas,
    pair_ids=pair_ids
  )
}

maximum_pau_sum_error <- function(pau, gene_ids) {
  errors <- numeric(0)
  for (gene in unique(gene_ids)) {
    values <- colSums(pau[gene_ids == gene, , drop=FALSE], na.rm=TRUE)
    available <- colSums(!is.na(pau[gene_ids == gene, , drop=FALSE])) > 0
    errors <- c(errors, abs(values[available] - 1))
  }
  if (length(errors)) max(errors) else NA_real_
}

maximum_effect_sum_error <- function(delta, gene_ids) {
  errors <- vapply(unique(gene_ids), function(gene) {
    values <- delta[gene_ids == gene]
    if (!any(!is.na(values))) return(NA_real_)
    abs(sum(values, na.rm=TRUE))
  }, numeric(1))
  errors <- errors[is.finite(errors)]
  if (length(errors)) max(errors) else NA_real_
}

finite_max <- function(values) {
  values <- values[is.finite(values)]
  if (length(values)) max(values) else NA_real_
}

scalar_text <- function(value, default="") {
  if (length(value) == 0 || is.na(value) || !nzchar(as.character(value))) default else as.character(value)
}

if ("--self-test" %in% args) {
  fixture <- matrix(
    c(75, 7500, 0, 50, 5000, 0, 25, 2500, 0, 50, 5000, 0),
    nrow=2, byrow=TRUE,
    dimnames=list(c("p1", "p2"), c("d1", "d2", "d0", "n1", "n2", "n0"))
  )
  genes <- c("g1", "g1")
  pau <- within_gene_pau(fixture, genes)
  stopifnot(all.equal(unname(colSums(pau[, c("d1", "d2", "n1", "n2")])), rep(1, 4)))
  stopifnot(all(is.na(pau[, c("d0", "n0")])))
  samples <- data.frame(
    sample_id=colnames(fixture),
    condition=c("D", "D", "D", "N", "N", "N"),
    subject=c("S1", "S2", "S0", "S1", "S2", "S0"),
    stringsAsFactors=FALSE,
    row.names=colnames(fixture)
  )
  paired <- paired_effect(pau, samples, "D", "N", "subject")
  stopifnot(all.equal(unname(paired$delta), c(-0.25, 0.25), tolerance=1e-12))
  unpaired <- unpaired_effect(pau, samples, "D", "N")
  stopifnot(all.equal(unname(unpaired$delta), c(-0.25, 0.25), tolerance=1e-12))
  stopifnot(maximum_pau_sum_error(pau, genes) < 1e-12)
  stopifnot(maximum_effect_sum_error(paired$delta, genes) < 1e-12)
  translated <- sum(paired$delta * c(100, 200))
  stopifnot(abs(translated - sum(paired$delta * c(10100, 10200))) < 1e-10)
  reverse_shift <- sum(paired$delta * -c(100, 200))
  stopifnot(sign(translated) == -sign(reverse_shift))
  permuted <- c(2, 1)
  stopifnot(all.equal(
    paired$delta[permuted],
    paired_effect(pau[permuted, , drop=FALSE], samples, "D", "N", "subject")$delta,
    tolerance=1e-12
  ))
  cat("APA-A2 raw-count PAU, paired-effect, strand, and invariance self-test: PASS\n")
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
min_abs_delta <- as.numeric(get_arg("--min-abs-delta-pau", "0.10"))
default_design_text <- get_arg("--design", "~ condition")
contrast_id <- get_arg("--contrast-id", "")
index_file <- get_arg("--index-file", file.path(outdir, "result_index.tsv"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

required_catalog <- c("pas_id", "gene_id", "strand", "summit_start", "assignment_status")
if (!all(required_catalog %in% colnames(catalog))) stop("Active-PAS catalog contract is incomplete")
if (anyDuplicated(catalog$pas_id)) stop("Active-PAS catalog pas_id values must be unique")
rownames(counts_df) <- counts_df$pas_id
sample_columns <- intersect(samples_all$sample_id, setdiff(colnames(counts_df), "pas_id"))
count_matrix <- as.matrix(counts_df[, sample_columns, drop=FALSE])
storage.mode(count_matrix) <- "integer"
if (nzchar(contrast_id)) {
  contrasts <- contrasts[contrasts$contrast_id == contrast_id, , drop=FALSE]
  if (nrow(contrasts) != 1) stop(paste("Expected exactly one contrast:", contrast_id))
}

index <- list()
for (i in seq_len(nrow(contrasts))) {
  con <- contrasts[i, ]
  keep <- samples_all$sample_id %in% sample_columns &
    samples_all$condition %in% c(con$denominator, con$numerator)
  sample_data <- samples_all[keep, , drop=FALSE]
  rownames(sample_data) <- sample_data$sample_id
  sample_data$condition <- relevel(factor(sample_data$condition), ref=con$denominator)
  design_text <- if (
    "resolved_design" %in% colnames(contrasts) && !is.na(con$resolved_design) && nzchar(con$resolved_design)
  ) con$resolved_design else default_design_text
  formulas <- dexseq_formulas(design_text)
  missing <- setdiff(formulas$variables, colnames(sample_data))
  if (length(missing)) stop("Design columns missing for ", con$contrast_id, ": ", paste(missing, collapse=", "))
  for (variable in formulas$variables) sample_data[[variable]] <- factor(sample_data[[variable]])
  design_matrix <- model.matrix(as.formula(design_text), data=sample_data)
  if (qr(design_matrix)$rank < ncol(design_matrix)) stop("Pair-specific design is not full rank: ", con$contrast_id)

  eligible <- catalog$assignment_status == "unique" & !is.na(catalog$gene_id) & nzchar(catalog$gene_id)
  gene_counts_first <- table(catalog$gene_id[eligible])
  eligible <- eligible & catalog$gene_id %in% names(gene_counts_first[gene_counts_first >= 2])
  eligible <- eligible & rowSums(count_matrix[catalog$pas_id, rownames(sample_data), drop=FALSE]) >= min_count
  gene_counts_second <- table(catalog$gene_id[eligible])
  eligible <- eligible & catalog$gene_id %in% names(gene_counts_second[gene_counts_second >= 2])
  selected <- as.character(catalog$pas_id[eligible])
  if (length(selected) < 2) stop("No testable APA-A2 sites for ", con$contrast_id)
  mat <- count_matrix[selected, rownames(sample_data), drop=FALSE]
  group_ids <- as.character(catalog$gene_id[match(selected, catalog$pas_id)])

  dxd <- DEXSeqDataSet(
    mat, sampleData=sample_data, design=formulas$full,
    featureID=selected, groupID=group_ids
  )
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd, reducedModel=formulas$reduced)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition")
  dexseq_result <- DEXSeqResults(dxd)
  gene_qvalues <- perGeneQValue(dexseq_result)
  result <- tabular_dexseq_result(as.data.frame(dexseq_result))
  result$pas_id <- as.character(result$featureID)
  catalog_rows <- catalog_rows_for_pas(catalog, result$pas_id)

  pau <- within_gene_pau(mat, group_ids)
  pairing_column <- if ("pairing_column" %in% colnames(contrasts)) scalar_text(con$pairing_column) else ""
  design_mode <- if ("design_mode" %in% colnames(contrasts)) scalar_text(con$design_mode, "unpaired") else "unpaired"
  effects <- if (design_mode == "paired") {
    if (!nzchar(pairing_column)) stop("Paired APA-A2 contrast lacks pairing_column: ", con$contrast_id)
    paired_effect(pau, sample_data, con$denominator, con$numerator, pairing_column)
  } else unpaired_effect(pau, sample_data, con$denominator, con$numerator)

  effect_rows <- match(result$pas_id, selected)
  result$gene_id <- as.character(catalog$gene_id[catalog_rows])
  result$gene_padj <- as.numeric(gene_qvalues[match(result$gene_id, names(gene_qvalues))])
  result$summit_start <- as.numeric(catalog$summit_start[catalog_rows])
  result$strand <- as.character(catalog$strand[catalog_rows])
  optional_catalog <- intersect(
    c("chrom", "summit_end", "start", "end", "feature_class", "assignment_status", "method", "interpretation"),
    colnames(catalog)
  )
  for (column in optional_catalog) result[[column]] <- catalog[[column]][catalog_rows]
  normalized <- counts(dxd, normalized=TRUE)
  rownames(normalized) <- selected
  result$mean_normalized_count <- rowMeans(normalized)[effect_rows]
  result$mean_raw_count <- rowMeans(mat)[effect_rows]
  result$PAU_denominator <- effects$denominator[effect_rows]
  result$PAU_numerator <- effects$numerator[effect_rows]
  result$delta_PAU <- effects$delta[effect_rows]
  result$significant_site <- !is.na(result$padj) & result$padj <= fdr
  result$primary_site <- result$significant_site & !is.na(result$delta_PAU) & abs(result$delta_PAU) >= min_abs_delta
  leading <- c(
    "pas_id", "gene_id", "chrom", "summit_start", "summit_end", "strand", "feature_class",
    "assignment_status", "gene_padj", "PAU_denominator", "PAU_numerator", "delta_PAU",
    "significant_site", "primary_site", "mean_raw_count", "mean_normalized_count"
  )
  leading <- intersect(leading, colnames(result))
  result <- result[, c(leading, setdiff(colnames(result), leading)), drop=FALSE]
  site_target <- file.path(outdir, paste0(con$contrast_id, ".apa_a2_sites.tsv"))
  write.table(result, site_target, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  pair_target <- file.path(outdir, paste0(con$contrast_id, ".apa_a2_pair_deltas.tsv"))
  pair_rows <- data.frame(
    pas_id=character(0), gene_id=character(0), pair_id=character(0),
    PAU_denominator=numeric(0), PAU_numerator=numeric(0), delta_PAU=numeric(0),
    stringsAsFactors=FALSE
  )
  if (length(effects$pair_ids)) {
    pair_rows <- do.call(rbind, lapply(seq_along(effects$pair_ids), function(pair_index) {
      den_sample <- rownames(sample_data)[
        sample_data$condition == con$denominator &
          as.character(sample_data[[pairing_column]]) == effects$pair_ids[pair_index]
      ]
      num_sample <- rownames(sample_data)[
        sample_data$condition == con$numerator &
          as.character(sample_data[[pairing_column]]) == effects$pair_ids[pair_index]
      ]
      data.frame(
        pas_id=selected, gene_id=group_ids, pair_id=effects$pair_ids[pair_index],
        PAU_denominator=pau[, den_sample], PAU_numerator=pau[, num_sample],
        delta_PAU=effects$pair_deltas[, pair_index], stringsAsFactors=FALSE
      )
    }))
  }
  write.table(pair_rows, pair_target, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  gene_summary <- do.call(rbind, lapply(split(result, result$gene_id), function(rows) {
    delta <- as.numeric(rows$delta_PAU)
    transcript_coordinate <- ifelse(rows$strand == "+", rows$summit_start, -rows$summit_start)
    weighted_shift <- if (any(!is.na(delta))) sum(delta * transcript_coordinate, na.rm=TRUE) else NA_real_
    significant_sites <- sum(rows$significant_site, na.rm=TRUE)
    primary_sites <- sum(rows$primary_site, na.rm=TRUE)
    gene_padj <- as.numeric(rows$gene_padj[1])
    significant_gene <- !is.na(gene_padj) && gene_padj <= fdr
    primary_gene <- significant_gene && primary_sites > 0
    shift <- if (!primary_gene || is.na(weighted_shift)) "no_shift" else if (weighted_shift > 0) {
      "distal"
    } else if (weighted_shift < 0) "proximal" else "no_directional_change"
    ordered <- order(transcript_coordinate)
    score <- if (significant_gene) -log10(max(gene_padj, .Machine$double.xmin)) else 0
    if (shift == "proximal") score <- -score
    if (!shift %in% c("proximal", "distal")) score <- 0
    data.frame(
      gene_id=rows$gene_id[1], gene_padj=gene_padj, testable_sites=nrow(rows),
      significant_sites=significant_sites, primary_sites=primary_sites,
      max_abs_delta_PAU=finite_max(abs(delta)),
      weighted_transcript_position_shift_nt=weighted_shift,
      shift=shift, proximal_pas=rows$pas_id[ordered[1]], distal_pas=rows$pas_id[ordered[length(ordered)]],
      significant_gene=significant_gene, primary_gene=primary_gene,
      signed_shift_score=score, stringsAsFactors=FALSE
    )
  }))
  gene_target <- file.path(outdir, paste0(con$contrast_id, ".apa_a2_genes.tsv"))
  write.table(gene_summary, gene_target, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  shift_target <- file.path(outdir, paste0(con$contrast_id, ".apa_a2_shifts.tsv"))
  write.table(
    gene_summary[, c(
      "gene_id", "gene_padj", "shift", "weighted_transcript_position_shift_nt",
      "proximal_pas", "distal_pas", "primary_sites", "primary_gene"
    )],
    shift_target, sep="\t", quote=FALSE, row.names=FALSE, na="NA"
  )

  zero_total_gene_samples <- sum(vapply(unique(group_ids), function(gene) {
    sum(colSums(mat[group_ids == gene, , drop=FALSE]) == 0)
  }, integer(1)))
  audit <- data.frame(
    contrast_id=con$contrast_id, status="PASS", design_mode=design_mode,
    effect_method=if (design_mode == "paired") "equal_weight_within_pair_raw_count_PAU" else "condition_mean_raw_count_PAU",
    matched_pairs=length(effects$pair_ids), selected_sites=length(selected),
    selected_genes=length(unique(group_ids)), zero_total_gene_samples=zero_total_gene_samples,
    maximum_PAU_sum_error=maximum_pau_sum_error(pau, group_ids),
    maximum_gene_delta_sum_error=maximum_effect_sum_error(effects$delta, group_ids),
    unavailable_site_effects=sum(is.na(effects$delta)),
    fdr=fdr, min_abs_delta_PAU=min_abs_delta, stringsAsFactors=FALSE
  )
  audit_target <- file.path(outdir, paste0(con$contrast_id, ".apa_a2_audit.tsv"))
  write.table(audit, audit_target, sep="\t", quote=FALSE, row.names=FALSE, na="NA")

  index[[length(index) + 1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=site_target, shift_file=shift_target,
    gene_summary_file=gene_target, pair_delta_file=pair_target, audit_file=audit_target,
    tested_sites=nrow(result), significant_sites=sum(result$significant_site, na.rm=TRUE),
    primary_sites=sum(result$primary_site, na.rm=TRUE), tested_genes=nrow(gene_summary),
    significant_genes=sum(gene_summary$significant_gene), primary_genes=sum(gene_summary$primary_gene),
    design_mode=design_mode, resolved_design=design_text,
    paired=if ("paired" %in% colnames(contrasts)) con$paired else design_mode == "paired",
    n_pairs=if ("n_pairs" %in% colnames(contrasts)) con$n_pairs else length(effects$pair_ids),
    warning=if ("design_status" %in% colnames(contrasts)) con$design_status else "",
    check.names=FALSE, stringsAsFactors=FALSE
  )
}

write.table(do.call(rbind, index), index_file, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
