# Regenerate only the distance figures from a saved matrix; no DESeq2 refit.
args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2) stop("Usage: Rscript replot_vst_distances.R DISTANCES.tsv OUTPUT_PREFIX")
suppressPackageStartupMessages(library(ggplot2))
script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
script_path <- normalizePath(sub("^--file=","",script_arg[[1]]),mustWork=TRUE)
plot_source <- file.path(dirname(script_path),"R","deseq2_gex.R")
expressions <- parse(plot_source)
selected <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) &&
  identical(x[[2]],as.name("sample_distance_plot")),as.list(expressions))
if (length(selected)!=1) stop("Could not locate the installed distance-plot function")
eval(selected[[1]])
distances <- as.matrix(read.delim(args[1],row.names=1,check.names=FALSE))
if (!is.numeric(distances) || nrow(distances)<2 || nrow(distances)!=ncol(distances) ||
    anyDuplicated(rownames(distances)) || anyDuplicated(colnames(distances)) ||
    !setequal(rownames(distances),colnames(distances)) || any(!is.finite(distances)) || any(distances<0)) {
  stop("Expected a finite, nonnegative square distance matrix with matching sample IDs")
}
distances <- distances[,rownames(distances),drop=FALSE]
if (!isSymmetric(distances,tol=1e-8) || any(abs(diag(distances))>1e-8)) stop("Distances must be symmetric with a zero diagonal")
prefix <- args[2]
dir.create(dirname(prefix),recursive=TRUE,showWarnings=FALSE)
plot <- sample_distance_plot(distances)
width <- max(8,3+0.28*ncol(distances)); height <- max(7,2+0.28*nrow(distances))
ggsave(paste0(prefix,".png"),plot,width=width,height=height,dpi=150)
ggsave(paste0(prefix,".pdf"),plot,width=width,height=height)
order <- levels(plot$data$sample_x)
write.table(data.frame(position=seq_along(order),sample_id=order),paste0(prefix,".sample_order.tsv"),
            sep="\t",quote=FALSE,row.names=FALSE)
cat("Regenerated clustered PNG/PDF; counts, distances and DE results were not recomputed.\n")
