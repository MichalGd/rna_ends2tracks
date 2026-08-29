args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args); if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name)); args[[pos + 1]]
}
suppressPackageStartupMessages({library(DRIMSeq); library(stageR)})

stagewise_with_na_policy <- function(gene_result, feature_result, alpha) {
  gene_ids <- as.character(gene_result$gene_id)
  feature_ids <- as.character(feature_result$feature_id)
  feature_gene_ids <- as.character(feature_result$gene_id)
  gene_pvalues <- as.numeric(gene_result$pvalue)
  feature_pvalues <- as.numeric(feature_result$pvalue)

  if (anyDuplicated(gene_ids)) stop("DRIMSeq returned duplicate gene IDs")
  if (anyDuplicated(feature_ids)) stop("DRIMSeq returned duplicate feature IDs")
  invalid_gene <- !is.na(gene_pvalues) & (!is.finite(gene_pvalues) | gene_pvalues < 0 | gene_pvalues > 1)
  invalid_feature <- !is.na(feature_pvalues) & (!is.finite(feature_pvalues) | feature_pvalues < 0 | feature_pvalues > 1)
  if (any(invalid_gene) || any(invalid_feature)) {
    stop("DRIMSeq returned finite p-values outside [0,1]")
  }

  screen_testable <- !is.na(gene_pvalues) & is.finite(gene_pvalues)
  confirmation_testable <- !is.na(feature_pvalues) & is.finite(feature_pvalues)
  finite_screen_genes <- gene_ids[screen_testable]
  testable_per_gene <- table(feature_gene_ids[
    confirmation_testable & feature_gene_ids %in% finite_screen_genes
  ])
  eligible_genes <- names(testable_per_gene)[testable_per_gene >= 2]
  eligible_features <- confirmation_testable & feature_gene_ids %in% eligible_genes
  eligible_screen <- screen_testable & gene_ids %in% eligible_genes

  if (!any(eligible_screen) || !any(eligible_features)) {
    stop("No genes with a finite screening p-value and at least two finite PAS confirmation p-values")
  }

  p_screen <- gene_pvalues[eligible_screen]
  names(p_screen) <- gene_ids[eligible_screen]
  p_confirmation <- matrix(
    feature_pvalues[eligible_features], ncol=1,
    dimnames=list(feature_ids[eligible_features], "site")
  )
  tx2gene <- data.frame(
    txID=feature_ids[eligible_features],
    geneID=feature_gene_ids[eligible_features],
    stringsAsFactors=FALSE
  )
  staged <- stageRTx(
    pScreen=p_screen, pConfirmation=p_confirmation,
    pScreenAdjusted=FALSE, tx2gene=tx2gene
  )
  staged <- stageWiseAdjustment(staged, method="dtu", alpha=alpha, allowNA=TRUE)
  adjusted <- as.data.frame(
    getAdjustedPValues(staged, order=FALSE, onlySignificantGenes=FALSE)
  )
  id_column <- intersect(c("txID", "feature_id"), colnames(adjusted))
  adjusted_ids <- if (length(id_column)) as.character(adjusted[[id_column[1]]]) else rownames(adjusted)
  value_column <- intersect(
    c("transcript", "feature", "site", "adjusted_pvalue", "padj"),
    colnames(adjusted)
  )
  if (!length(value_column)) stop("stageR returned no feature-level adjusted-p-value column")
  adjusted_locations <- match(adjusted_ids, feature_ids)
  if (anyNA(adjusted_locations)) stop("stageR returned an unknown feature ID")
  stage_adjusted <- rep(NA_real_, length(feature_ids))
  names(stage_adjusted) <- feature_ids
  stage_adjusted[adjusted_locations] <- as.numeric(adjusted[[value_column[1]]])

  fraction <- function(numerator, denominator) {
    if (denominator) numerator / denominator else NA_real_
  }
  screen_na <- sum(!screen_testable)
  confirmation_na <- sum(!confirmation_testable)
  excluded_lt2 <- sum(screen_testable & !gene_ids %in% eligible_genes)
  adjusted_na <- sum(is.na(stage_adjusted))
  audit <- data.frame(
    status=ifelse(screen_na + confirmation_na + excluded_lt2 > 0, "WARN_UNTESTABLE_PVALUES", "PASS"),
    screening_tests=length(gene_pvalues),
    screening_na=screen_na,
    screening_na_fraction=fraction(screen_na, length(gene_pvalues)),
    confirmation_tests=length(feature_pvalues),
    confirmation_na=confirmation_na,
    confirmation_na_fraction=fraction(confirmation_na, length(feature_pvalues)),
    excluded_genes_fewer_than_two_testable_sites=excluded_lt2,
    stageR_input_genes=length(p_screen),
    stageR_input_sites=nrow(p_confirmation),
    stageR_adjusted_na=adjusted_na,
    stageR_adjusted_na_fraction=fraction(adjusted_na, length(stage_adjusted)),
    na_policy="untestable hypotheses remain NA and cannot be significant",
    stringsAsFactors=FALSE
  )
  if (audit$status != "PASS") {
    warning(sprintf(
      paste0(
        "Retaining untestable hypotheses as NA: screening=%d/%d; ",
        "confirmation=%d/%d; genes_with_fewer_than_two_testable_sites=%d"
      ),
      screen_na, length(gene_pvalues), confirmation_na, length(feature_pvalues), excluded_lt2
    ))
  }
  list(adjusted=unname(stage_adjusted), audit=audit)
}

contrast_seed <- function(contrast_id) {
  codepoints <- utf8ToInt(as.character(contrast_id))
  if (!length(codepoints)) return(104729L)
  as.integer((sum(codepoints * seq_along(codepoints)) %% 2147483000) + 1)
}

fit_drimseq_with_numeric_retry <- function(d, design, condition_coef, contrast_id, multifactor) {
  seed <- contrast_seed(contrast_id)
  one_way <- !isTRUE(multifactor)

  fit_once <- function(add_uniform) {
    set.seed(seed)
    fitted <- dmPrecision(
      d, design=design, one_way=one_way, add_uniform=add_uniform
    )
    set.seed(seed)
    fitted <- dmFit(
      fitted, design=design, one_way=one_way, add_uniform=add_uniform
    )
    dmTest(fitted, coef=condition_coef, one_way=one_way)
  }

  primary_error <- NULL
  fitted <- tryCatch(
    fit_once(FALSE),
    error=function(error) {
      primary_error <<- conditionMessage(error)
      NULL
    }
  )
  policy <- "standard"

  if (is.null(fitted)) {
    recognized_numeric_failure <- grepl(
      "BiocParallel errors|non-finite value supplied by optim|NaNs produced|optimHess",
      primary_error
    )
    if (!isTRUE(multifactor) || !recognized_numeric_failure) {
      stop(primary_error, call.=FALSE)
    }
    warning(sprintf(
      paste0(
        "DRIMSeq numerical zero-pattern failure for multifactor contrast %s; ",
        "retrying deterministically with documented add_uniform=TRUE"
      ),
      contrast_id
    ))
    fitted <- fit_once(TRUE)
    policy <- "deterministic_add_uniform_retry"
  }

  list(
    fitted=fitted,
    audit=data.frame(
      contrast_id=contrast_id,
      status=ifelse(policy == "standard", "PASS", "WARN_NUMERIC_RETRY"),
      fit_policy=policy,
      multifactor=isTRUE(multifactor),
      one_way=one_way,
      random_seed=seed,
      add_uniform_used=policy != "standard",
      primary_error=ifelse(is.null(primary_error), "", gsub("[\r\n\t]+", " ", primary_error)),
      stringsAsFactors=FALSE
    )
  )
}

if ("--self-test" %in% args) {
  self_gene <- data.frame(
    gene_id=c("g1", "g2", "g3", "g4"),
    pvalue=c(1e-8, NA_real_, 0.2, 0.9),
    stringsAsFactors=FALSE
  )
  self_feature <- data.frame(
    feature_id=paste0("f", 1:8),
    gene_id=rep(c("g1", "g2", "g3", "g4"), each=2),
    pvalue=c(1e-6, 0.2, 0.4, 0.5, 0.3, NA_real_, 0.8, 0.9),
    stringsAsFactors=FALSE
  )
  self <- suppressWarnings(stagewise_with_na_policy(self_gene, self_feature, 0.05))
  stopifnot(
    length(self$adjusted) == nrow(self_feature),
    all(is.na(self$adjusted[3:6])),
    self$audit$screening_na == 1,
    self$audit$confirmation_na == 1,
    self$audit$excluded_genes_fewer_than_two_testable_sites == 1,
    self$audit$stageR_input_genes == 2,
    self$audit$stageR_input_sites == 4,
    identical(self$audit$status, "WARN_UNTESTABLE_PVALUES")
  )
  stopifnot(
    identical(contrast_seed("CTCF_IAA_vs_CTCF_control"), contrast_seed("CTCF_IAA_vs_CTCF_control")),
    contrast_seed("CTCF_IAA_vs_CTCF_control") != contrast_seed("RAD21_control_vs_CTCF_control")
  )
  cat("DRIMSeq/stageR NA-policy and deterministic-fit self-test: PASS\n")
  quit(save="no", status=0)
}

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
  multifactor <- length(setdiff(design_variables, "condition")) > 0
  fit <- fit_drimseq_with_numeric_retry(
    d, design, condition_coef, as.character(con$contrast_id), multifactor
  )
  d <- fit$fitted
  gene_result <- DRIMSeq::results(d, level="gene"); feature_result <- DRIMSeq::results(d, level="feature")
  gene_result$gene_padj <- p.adjust(gene_result$pvalue, method="BH")
  stagewise <- stagewise_with_na_policy(gene_result, feature_result, alpha)
  feature_result$stageR_adjusted <- stagewise$adjusted
  den <- rowMeans(mat[, samples$condition == con$denominator, drop=FALSE])
  num <- rowMeans(mat[, samples$condition == con$numerator, drop=FALSE])
  den_total <- ave(den, gene_ids, FUN=sum); num_total <- ave(num, gene_ids, FUN=sum)
  feature_result$PAU_denominator <- den[match(feature_result$feature_id, rownames(mat))] / pmax(den_total[match(feature_result$feature_id, rownames(mat))], 1e-12)
  feature_result$PAU_numerator <- num[match(feature_result$feature_id, rownames(mat))] / pmax(num_total[match(feature_result$feature_id, rownames(mat))], 1e-12)
  feature_result$delta_PAU <- feature_result$PAU_numerator - feature_result$PAU_denominator
  positions <- catalog$start[match(feature_result$feature_id, catalog$pas_id)]
  result_genes <- feature_result$gene_id
  feature_result$weighted_genomic_position_shift_nt <- ave(feature_result$PAU_numerator * positions, result_genes, FUN=sum) - ave(feature_result$PAU_denominator * positions, result_genes, FUN=sum)
  feature_strands <- catalog$strand[match(feature_result$feature_id, catalog$pas_id)]
  feature_result$weighted_transcript_position_shift_nt <- feature_result$weighted_genomic_position_shift_nt * ifelse(feature_strands == "+", 1, -1)
  finite_max <- function(values) {
    values <- values[is.finite(values)]
    if (length(values)) max(values) else NA_real_
  }
  gene_summary <- do.call(rbind, lapply(split(feature_result, feature_result$gene_id), function(rows) {
    gene_id <- as.character(rows$gene_id[1])
    gene_padj <- gene_result$gene_padj[match(gene_id, gene_result$gene_id)]
    shift_nt <- rows$weighted_transcript_position_shift_nt[1]
    data.frame(
      gene_id=gene_id, gene_padj=gene_padj, testable_sites=nrow(rows),
      confirmed_sites=sum(!is.na(rows$stageR_adjusted) & rows$stageR_adjusted < alpha),
      max_abs_delta_PAU=finite_max(abs(rows$delta_PAU)),
      weighted_transcript_position_shift_nt=shift_nt,
      shift=ifelse(is.na(shift_nt) || shift_nt == 0, "no_directional_change", ifelse(shift_nt > 0, "distal", "proximal")),
      signed_shift_score=ifelse(is.na(gene_padj), NA_real_, sign(shift_nt) * -log10(max(gene_padj, .Machine$double.xmin))),
      stringsAsFactors=FALSE
    )
  }))
  target <- file.path(outdir, paste0(con$contrast_id, ".drimseq_stager.tsv"))
  gene_target <- file.path(outdir, paste0(con$contrast_id, ".gene_apa_summary.tsv"))
  na_audit_target <- file.path(outdir, paste0(con$contrast_id, ".na_audit.tsv"))
  fit_audit_target <- file.path(outdir, paste0(con$contrast_id, ".fit_audit.tsv"))
  write.table(feature_result, target, sep="\t", quote=FALSE, row.names=FALSE)
  gene_screen_target <- file.path(outdir, paste0(con$contrast_id, ".gene_screen.tsv"))
  write.table(gene_result, gene_screen_target, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(gene_summary, gene_target, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(stagewise$audit, na_audit_target, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(fit$audit, fit_audit_target, sep="\t", quote=FALSE, row.names=FALSE)
  index[[length(index)+1]] <- data.frame(
    contrast_id=con$contrast_id, result_file=target, tested_sites=nrow(feature_result),
    gene_screen_file=gene_screen_target, gene_summary_file=gene_target,
    na_audit_file=na_audit_target,
    fit_audit_file=fit_audit_target,
    fit_policy=fit$audit$fit_policy,
    screening_na=stagewise$audit$screening_na,
    confirmation_na=stagewise$audit$confirmation_na,
    na_audit_status=stagewise$audit$status,
    confirmed_sites=sum(!is.na(feature_result$stageR_adjusted) & feature_result$stageR_adjusted < alpha),
    design_mode=if ("design_mode" %in% colnames(contrasts)) con$design_mode else "legacy",
    resolved_design=design_text,
    paired=if ("paired" %in% colnames(contrasts)) con$paired else FALSE,
    n_pairs=if ("n_pairs" %in% colnames(contrasts)) con$n_pairs else 0,
    warning=con$design_status, check.names=FALSE
  )
}
write.table(do.call(rbind, index), index_file, sep="\t", quote=FALSE, row.names=FALSE)
