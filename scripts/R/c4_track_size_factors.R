args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(name) {
  pos <- match(name, args)
  if (is.na(pos) || pos == length(args)) stop(paste("Required argument:", name))
  args[[pos + 1]]
}
suppressPackageStartupMessages(library(DESeq2))
raw <- read.delim(get_arg("--counts"), check.names=FALSE, stringsAsFactors=FALSE)
samples <- read.delim(get_arg("--samples"), check.names=FALSE, stringsAsFactors=FALSE)
output <- get_arg("--output")
if (!"gene_id" %in% colnames(raw)) stop("C4 matrix requires gene_id")
sample_ids <- intersect(samples$sample_id, setdiff(colnames(raw), "gene_id"))
matrix <- as.matrix(raw[, sample_ids, drop=FALSE])
rownames(matrix) <- raw$gene_id
storage.mode(matrix) <- "integer"
metadata <- samples[match(sample_ids, samples$sample_id), , drop=FALSE]
rownames(metadata) <- metadata$sample_id
dds <- DESeqDataSetFromMatrix(matrix, metadata, ~ 1)
dds <- estimateSizeFactors(dds, type="poscounts")
size_factor <- sizeFactors(dds)
column_sum <- colSums(matrix)
positive <- column_sum[column_sum > 0]
geometric_mean <- if (length(positive)) exp(mean(log(positive))) else stop("All C4 libraries are zero")
table <- data.frame(
  sample_id=names(size_factor), size_factor=as.numeric(size_factor), C4_column_sum=as.numeric(column_sum),
  cohort_geometric_mean_column_sum=geometric_mean, deseq2_scale=1/as.numeric(size_factor),
  robust_effective_library=as.numeric(size_factor)*geometric_mean,
  robust_cpm_scale=1000000/(as.numeric(size_factor)*geometric_mean), estimator="DESeq2_poscounts",
  stringsAsFactors=FALSE
)
expected <- DESeq2::fpm(dds, robust=TRUE)
observed <- sweep(counts(dds), 2, table$robust_effective_library, "/") * 1000000
tolerance <- 1e-8 * max(1, max(abs(expected)))
if (max(abs(expected-observed)) > tolerance) stop("Robust CPM formula differs from DESeq2::fpm(..., robust=TRUE)")
dir.create(dirname(output), recursive=TRUE, showWarnings=FALSE)
write.table(table, output, sep="\t", quote=FALSE, row.names=FALSE)
