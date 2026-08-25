set.seed(20260825)
root <- normalizePath(file.path(tempdir(), "rna_ends2tracks_pairing_smoke"), mustWork=FALSE)
dir.create(root, recursive=TRUE, showWarnings=FALSE)
outdir <- file.path(root, "results")

samples <- data.frame(
  sample_id=c("A1", "A2", "A3", "B1", "B2", "B3", "C1", "C2", "C3"),
  condition=c(rep("A", 3), rep("B", 3), rep("C", 3)),
  subject=c("S1", "S2", "S3", "S1", "S2", "S3", "C1", "C2", "C3"),
  stringsAsFactors=FALSE
)
samples_path <- file.path(root, "samples.tsv")
write.table(samples, samples_path, sep="\t", quote=FALSE, row.names=FALSE)

contrasts <- data.frame(
  contrast_id=c("B_vs_A", "C_vs_A"),
  factor="condition",
  numerator=c("B", "C"),
  denominator="A",
  n_num=3,
  n_den=3,
  paired=c(TRUE, FALSE),
  n_pairs=c(3, 0),
  pairing_status=c("complete", "disjoint_subjects"),
  design_status="valid",
  design_mode=c("paired", "unpaired"),
  resolved_design=c("~ subject + condition", "~ condition"),
  pairing_column=c("subject", ""),
  stringsAsFactors=FALSE
)
contrasts_path <- file.path(root, "contrasts.tsv")
write.table(contrasts, contrasts_path, sep="\t", quote=FALSE, row.names=FALSE)

gene_count <- 1200L
sample_count <- nrow(samples)
matrix <- matrix(rnbinom(gene_count * sample_count, mu=100, size=20), nrow=gene_count)
matrix[1:100, samples$condition == "B"] <- matrix[1:100, samples$condition == "B"] + 100L
matrix[101:200, samples$condition == "C"] <- matrix[101:200, samples$condition == "C"] + 80L
colnames(matrix) <- samples$sample_id
counts <- data.frame(
  Geneid=paste0("gene", seq_len(gene_count)), Chr="chr1", Start=seq_len(gene_count),
  End=seq_len(gene_count), Strand="+", Length=100L, matrix, check.names=FALSE
)
counts_path <- file.path(root, "featureCounts.tsv")
write.table(counts, counts_path, sep="\t", quote=FALSE, row.names=FALSE)

status <- system2("Rscript", c(
  "scripts/R/deseq2_all_pairs.R",
  "--counts", counts_path,
  "--samples", samples_path,
  "--contrasts", contrasts_path,
  "--design", shQuote("~ condition"),
  "--outdir", outdir,
  "--fdr", "0.05"
))
stopifnot(status == 0L)
index <- read.delim(file.path(outdir, "result_index.tsv"), check.names=FALSE, stringsAsFactors=FALSE)
stopifnot(nrow(index) == 2L)
stopifnot(index$design_mode == c("paired", "unpaired"))
stopifnot(index$resolved_design == c("~ subject + condition", "~ condition"))
stopifnot(index$n_pairs == c(3L, 0L))
stopifnot(all(file.exists(index$result_file)))
cat("DESeq2 mixed pairing smoke: PASS\n")
