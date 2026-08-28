args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name, default=NULL) {
  pos <- match(name, args)
  if (is.na(pos)) return(default)
  if (pos == length(args)) stop(paste("Missing value for", name))
  args[[pos + 1]]
}
as_flag <- function(value) tolower(value) == "true"
suppressPackageStartupMessages({library(msigdbr); library(fgsea); library(ggplot2)})
if ("--self-test" %in% args) {
  fetch <- function(species, collection, subcollection=NULL) {
    formal_names <- names(formals(msigdbr::msigdbr))
    call <- list(species=species)
    if ("db_species" %in% formal_names) call$db_species <- "HS"
    if ("collection" %in% formal_names) call$collection <- collection else call$category <- collection
    if (!is.null(subcollection)) {
      if ("subcollection" %in% formal_names) call$subcollection <- subcollection else call$subcategory <- subcollection
    }
    suppressWarnings(do.call(msigdbr::msigdbr, call))
  }
  hallmark_human <- fetch("Homo sapiens", "H")
  hallmark_mouse <- fetch("Mus musculus", "H")
  go_bp <- fetch("Homo sapiens", "C5", "GO:BP")
  reactome <- fetch("Homo sapiens", "C2", "CP:REACTOME")
  stopifnot(nrow(hallmark_human) > 0, nrow(hallmark_mouse) > 0, nrow(go_bp) > 0, nrow(reactome) > 0)
  gene_column <- if ("ensembl_gene" %in% colnames(hallmark_human)) "ensembl_gene" else "db_ensembl_gene"
  pathways <- split(hallmark_human[[gene_column]], hallmark_human$gs_id)
  gene_ids <- unique(unlist(pathways)); gene_ids <- gene_ids[nzchar(gene_ids)]
  stats <- setNames(seq(-2, 2, length.out=length(gene_ids)), gene_ids)
  result <- suppressWarnings(fgsea::fgseaMultilevel(pathways=head(pathways, 3), stats=stats, minSize=5, maxSize=500))
  stopifnot(is.data.frame(as.data.frame(result)))
  invisible(ggplot2::ggplot(data.frame(x=1, y=1), ggplot2::aes(x, y)) + ggplot2::geom_point())
  message("Gene-set database and fgsea self-test: PASS")
  quit(save="no", status=0)
}
required <- c("--input", "--outdir", "--analysis-type", "--species", "--genome", "--contrast-id")
for (item in required) if (is.null(get_arg(item))) stop(paste("Required argument:", item))

input <- get_arg("--input")
outdir <- get_arg("--outdir")
analysis_type <- get_arg("--analysis-type")
species <- get_arg("--species")
target_species <- if (tolower(species) == "mouse") "Mus musculus" else "Homo sapiens"
genome <- get_arg("--genome")
contrast_id <- get_arg("--contrast-id")
run_ora <- as_flag(get_arg("--ora", "true"))
run_gsea <- as_flag(get_arg("--gsea", "true"))
use_go <- as_flag(get_arg("--go", "true"))
use_reactome <- as_flag(get_arg("--reactome", "true"))
use_hallmarks <- as_flag(get_arg("--hallmarks", "true"))
padj_threshold <- as.numeric(get_arg("--padj", "0.05"))
min_size <- as.integer(get_arg("--min-size", "10"))
max_size <- as.integer(get_arg("--max-size", "500"))
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

genes <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
required_columns <- c("gene_id", "background", "foreground", "direction", "rank_score")
if (!all(required_columns %in% colnames(genes))) stop("Prepared enrichment table has an invalid schema")
genes$gene_id <- sub("[.][0-9]+$", "", as.character(genes$gene_id))
genes <- genes[!is.na(genes$gene_id) & nzchar(genes$gene_id), , drop=FALSE]
genes <- genes[!duplicated(genes$gene_id), , drop=FALSE]

msig_call <- function(collection, subcollection=NULL, database) {
  formal_names <- names(formals(msigdbr::msigdbr))
  call <- list(species=target_species)
  if ("db_species" %in% formal_names) call$db_species <- "HS"
  if ("collection" %in% formal_names) call$collection <- collection else call$category <- collection
  if (!is.null(subcollection)) {
    if ("subcollection" %in% formal_names) call$subcollection <- subcollection else call$subcategory <- subcollection
  }
  frame <- suppressWarnings(do.call(msigdbr::msigdbr, call))
  if (!nrow(frame)) return(data.frame())
  gene_column <- if ("ensembl_gene" %in% colnames(frame)) "ensembl_gene" else "db_ensembl_gene"
  data.frame(
    database=database, term_id=as.character(frame$gs_id), term_name=as.character(frame$gs_name),
    gene_id=sub("[.][0-9]+$", "", as.character(frame[[gene_column]])),
    db_version=if ("db_version" %in% colnames(frame)) as.character(frame$db_version) else NA_character_,
    stringsAsFactors=FALSE
  )
}

sets <- list()
if (use_go) {
  sets[[length(sets) + 1]] <- msig_call("C5", "GO:BP", "GO_BP")
  sets[[length(sets) + 1]] <- msig_call("C5", "GO:MF", "GO_MF")
  sets[[length(sets) + 1]] <- msig_call("C5", "GO:CC", "GO_CC")
}
if (use_reactome) sets[[length(sets) + 1]] <- msig_call("C2", "CP:REACTOME", "REACTOME")
if (use_hallmarks) sets[[length(sets) + 1]] <- msig_call("H", NULL, "HALLMARK")
sets <- sets[vapply(sets, nrow, integer(1)) > 0]
term2gene <- if (length(sets)) unique(do.call(rbind, sets)) else data.frame(
  database=character(), term_id=character(), term_name=character(), gene_id=character(), db_version=character()
)
term2gene <- term2gene[!is.na(term2gene$gene_id) & nzchar(term2gene$gene_id), , drop=FALSE]

ora_columns <- c("analysis_type", "contrast_id", "genome", "query", "database", "term_id", "term_name",
                 "set_size", "background_size", "foreground_size", "overlap_count", "pvalue", "padj", "genes")
ora_rows <- list()
background <- unique(genes$gene_id[as.logical(genes$background)])
direction_values <- strsplit(as.character(genes$direction), ";", fixed=TRUE)
queries <- sort(unique(unlist(direction_values[as.logical(genes$foreground)])))
queries <- queries[nzchar(queries)]
has_query <- function(query) vapply(direction_values, function(values) query %in% values, logical(1))
if (run_ora && length(background) && nrow(term2gene)) {
  for (query in queries) {
    foreground <- intersect(unique(genes$gene_id[as.logical(genes$foreground) & has_query(query)]), background)
    for (database in unique(term2gene$database)) {
      db <- term2gene[term2gene$database == database & term2gene$gene_id %in% background, , drop=FALSE]
      for (term_id in unique(db$term_id)) {
        term <- db[db$term_id == term_id, , drop=FALSE]
        members <- unique(term$gene_id); overlap <- intersect(foreground, members)
        if (length(members) < min_size || length(members) > max_size || !length(foreground)) next
        pvalue <- phyper(length(overlap) - 1, length(members), length(background) - length(members),
                         length(foreground), lower.tail=FALSE)
        ora_rows[[length(ora_rows) + 1]] <- data.frame(
          analysis_type=analysis_type, contrast_id=contrast_id, genome=genome, query=query,
          database=database, term_id=term_id, term_name=term$term_name[1], set_size=length(members),
          background_size=length(background), foreground_size=length(foreground),
          overlap_count=length(overlap), pvalue=pvalue, padj=NA_real_,
          genes=paste(sort(overlap), collapse=";"), stringsAsFactors=FALSE
        )
      }
    }
  }
}
ora <- if (length(ora_rows)) do.call(rbind, ora_rows) else as.data.frame(setNames(replicate(length(ora_columns), logical(0), simplify=FALSE), ora_columns))
if (nrow(ora)) ora$padj <- ave(ora$pvalue, interaction(ora$query, ora$database), FUN=function(x) p.adjust(x, "BH"))
write.table(ora, file.path(outdir, "ora.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

gsea_columns <- c("analysis_type", "contrast_id", "genome", "database", "term_id", "term_name",
                  "set_size", "ES", "NES", "pvalue", "padj", "leading_edge")
gsea_rows <- list()
ranked <- genes[as.logical(genes$background) & is.finite(genes$rank_score), c("gene_id", "rank_score"), drop=FALSE]
ranked <- ranked[order(-abs(ranked$rank_score)), , drop=FALSE]
ranked <- ranked[!duplicated(ranked$gene_id), , drop=FALSE]
stats <- sort(setNames(ranked$rank_score, ranked$gene_id), decreasing=TRUE)
if (run_gsea && length(stats) >= min_size && length(unique(stats)) > 1 && nrow(term2gene)) {
  set.seed(1)
  for (database in unique(term2gene$database)) {
    db <- term2gene[term2gene$database == database, , drop=FALSE]
    pathways <- split(db$gene_id, db$term_id)
    pathways <- lapply(pathways, unique)
    fg <- suppressWarnings(fgsea::fgseaMultilevel(
      pathways=pathways, stats=stats, minSize=min_size, maxSize=max_size, eps=0
    ))
    if (!nrow(fg)) next
    names_by_id <- setNames(db$term_name, db$term_id)
    frame <- as.data.frame(fg)
    frame$database <- database
    frame$term_id <- frame$pathway
    frame$term_name <- unname(names_by_id[frame$pathway])
    frame$leading_edge <- vapply(frame$leadingEdge, function(x) paste(x, collapse=";"), character(1))
    gsea_rows[[length(gsea_rows) + 1]] <- data.frame(
      analysis_type=analysis_type, contrast_id=contrast_id, genome=genome,
      database=frame$database, term_id=frame$term_id, term_name=frame$term_name,
      set_size=frame$size, ES=frame$ES, NES=frame$NES, pvalue=frame$pval,
      padj=frame$padj, leading_edge=frame$leading_edge, stringsAsFactors=FALSE
    )
  }
}
gsea <- if (length(gsea_rows)) do.call(rbind, gsea_rows) else as.data.frame(setNames(replicate(length(gsea_columns), logical(0), simplify=FALSE), gsea_columns))
write.table(gsea, file.path(outdir, "gsea.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

mapping <- data.frame(
  metric=c("input_genes", "background_genes", "mapped_background_genes", "mapping_fraction", "foreground_genes"),
  value=c(nrow(genes), length(background), length(intersect(background, term2gene$gene_id)),
          ifelse(length(background), length(intersect(background, term2gene$gene_id)) / length(background), NA_real_),
          sum(as.logical(genes$foreground))), stringsAsFactors=FALSE
)
write.table(mapping, file.path(outdir, "mapping_audit.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

plot_rows <- data.frame()
if (nrow(ora)) {
  selected <- ora[!is.na(ora$padj) & ora$padj <= padj_threshold, , drop=FALSE]
  if (nrow(selected)) {
    selected$score <- -log10(pmax(selected$padj, .Machine$double.xmin))
    selected$source <- paste("ORA", selected$query, selected$database, sep=":")
    plot_rows <- rbind(plot_rows, selected[, c("term_name", "score", "source")])
  }
}
if (nrow(gsea)) {
  selected <- gsea[!is.na(gsea$padj) & gsea$padj <= padj_threshold, , drop=FALSE]
  if (nrow(selected)) {
    selected$score <- abs(selected$NES)
    selected$source <- paste("GSEA", selected$database, sep=":")
    plot_rows <- rbind(plot_rows, selected[, c("term_name", "score", "source")])
  }
}
if (nrow(plot_rows)) {
  plot_rows <- plot_rows[order(-plot_rows$score), , drop=FALSE]
  plot_rows <- head(plot_rows, 25)
  plot_rows$term_name <- factor(plot_rows$term_name, levels=rev(unique(plot_rows$term_name)))
  enrichment_plot <- ggplot(plot_rows, aes(score, term_name, color=source)) + geom_point(size=2.5) +
    labs(x="Enrichment score", y=NULL, title=paste(analysis_type, contrast_id)) + theme_bw()
} else {
  enrichment_plot <- ggplot() + annotate("text", x=0, y=0, label="No enriched term passed the configured threshold") +
    xlim(-1, 1) + ylim(-1, 1) + theme_void() + labs(title=paste(analysis_type, contrast_id))
}
ggsave(file.path(outdir, "enrichment.pdf"), enrichment_plot, width=9, height=7)
ggsave(file.path(outdir, "enrichment.png"), enrichment_plot, width=9, height=7, dpi=150)

versions <- unique(na.omit(term2gene$db_version))
provenance <- data.frame(
  key=c("analysis_type", "contrast_id", "genome", "species", "msigdbr", "fgsea", "ggplot2",
        "msigdb_versions", "msigdb_source_species", "ora", "gsea", "padj", "min_size", "max_size"),
  value=c(analysis_type, contrast_id, genome, species, as.character(packageVersion("msigdbr")),
          as.character(packageVersion("fgsea")), as.character(packageVersion("ggplot2")),
          paste(versions, collapse=","), "Homo sapiens; mouse uses msigdbr ortholog mapping",
          run_ora, run_gsea, padj_threshold, min_size, max_size), stringsAsFactors=FALSE
)
write.table(provenance, file.path(outdir, "provenance.tsv"), sep="\t", quote=FALSE, row.names=FALSE)
