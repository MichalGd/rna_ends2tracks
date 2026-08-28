#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
value <- function(flag, default = NULL) {
  index <- match(flag, args)
  if (is.na(index)) return(default)
  if (index == length(args)) stop("Missing value after ", flag)
  args[[index + 1L]]
}

ends_file <- value("--ends")
fasta <- value("--fasta")
output <- value("--output")
audit_file <- value("--audit")
cluster_gap <- as.integer(value("--cluster-gap", "24"))
if (any(vapply(list(ends_file, fasta, output, audit_file), is.null, logical(1)))) {
  stop("Required: --ends --fasta --output --audit")
}

suppressPackageStartupMessages({
  library(PolyAseqTrap)
  library(Rsamtools)
  library(GenomicRanges)
  library(IRanges)
})

ends <- read.delim(ends_file, check.names = FALSE, stringsAsFactors = FALSE)
required <- c("chrom", "position", "strand", "count")
missing <- setdiff(required, names(ends))
if (length(missing)) stop("Exact-end table lacks columns: ", paste(missing, collapse = ", "))
if (!nrow(ends)) stop("Exact-end table is empty")
if (any(!ends$strand %in% c("+", "-")) || any(ends$count < 1L)) stop("Invalid exact-end rows")

# Python emits zero-based cleavage bases. PolyAseqTrap/GRanges uses one-based
# genomic positions, so add one here and subtract one again in the compact output.
points <- GenomicRanges::GRanges(
  seqnames = ends$chrom,
  ranges = IRanges::IRanges(start = as.integer(ends$position) + 1L, width = 1L),
  strand = ends$strand,
  score = as.numeric(ends$count)
)
clusters <- PolyAseqTrap::simpleCluster(points, max.gapwidth = cluster_gap)
# The pinned simpleCluster implementation assigns sum(extractList(...)) as a
# scalar and recycles that project-wide total to every cluster. Preserve its
# official strand-aware ranges and centers, but restore the intended per-PAC
# weighted count from each cluster's revmap members.
clusters$score <- vapply(
  seq_along(clusters$revmap),
  function(index) sum(points$score[clusters$revmap[[index]]]),
  numeric(1)
)
sites <- as.data.frame(clusters)

# check.repeat is part of the pinned PolyAseqTrap implementation and accepts a
# FaFile through Biostrings::getSeq. Edge sites that cannot supply ten bases are
# conservatively marked for DeepIP review in the Python adapter.
sites$coord <- as.integer(sites$center)
reference <- Rsamtools::FaFile(fasta)
Rsamtools::open.FaFile(reference)
on.exit(Rsamtools::close.FaFile(reference), add = TRUE)
repeat_flag <- tryCatch(
  PolyAseqTrap::check.repeat(sites, bsgenome = reference)$repeat_flag,
  error = function(condition) {
    warning("PolyAseqTrap repeat check deferred to adapter: ", conditionMessage(condition))
    rep(TRUE, nrow(sites))
  }
)

compact <- data.frame(
  chrom = as.character(sites$seqnames),
  position = as.integer(sites$center) - 1L,
  strand = as.character(sites$strand),
  count = as.integer(round(sites$score)),
  polyA_supported_count = NA_integer_,
  coordinate_level = "QuantSeq_REV_no_tail_weighted_PAC",
  repeat_detected = as.logical(repeat_flag),
  stringsAsFactors = FALSE
)
compact <- compact[order(compact$chrom, compact$position, compact$strand), , drop = FALSE]
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.table(compact, output, sep = "\t", quote = FALSE, row.names = FALSE)

audit <- data.frame(
  metric = c("eligible_endpoint_records", "distinct_exact_end_coordinates", "weighted_PACs"),
  value = c(sum(ends$count), nrow(ends), nrow(compact)),
  stringsAsFactors = FALSE
)
write.table(audit, audit_file, sep = "\t", quote = FALSE, row.names = FALSE)
