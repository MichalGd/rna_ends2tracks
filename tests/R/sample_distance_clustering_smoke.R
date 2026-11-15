# Focused plotting regression; evaluates only the actual plot function.
args <- commandArgs(trailingOnly=TRUE)
repo <- normalizePath(args[1],winslash="/",mustWork=TRUE)
output <- normalizePath(args[2],winslash="/",mustWork=FALSE)
dir.create(output,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages(library(ggplot2))
scripts <- file.path(repo,"scripts","R",c("deseq2_gex.R","deseq2_c4.R"))
load_plot <- function(path) {
  expressions <- parse(path)
  selected <- Filter(function(x) is.call(x) && identical(x[[1]],as.name("<-")) &&
    identical(x[[2]],as.name("sample_distance_plot")),as.list(expressions))
  stopifnot(length(selected)==1)
  environment <- new.env(parent=globalenv())
  eval(selected[[1]],environment)
  environment$sample_distance_plot
}
functions <- lapply(scripts,load_plot)
stopifnot(identical(body(functions[[1]]),body(functions[[2]])))
# Alphabetical order interleaves two well-separated biological groups.
vst <- matrix(c(0,0,100,100,0.2,0.3,100.3,100.2),nrow=2,
              dimnames=list(NULL,c("A_near","B_far","Y_near","Z_far")))
for (permutation in list(1:4,4:1)) {
  distances <- as.matrix(dist(t(vst[,permutation])))
  original <- distances
  for (function_index in seq_along(functions)) {
    plot <- functions[[function_index]](distances)
    order <- levels(plot$data$sample_x)
    stopifnot(identical(levels(plot$data$sample_y),rev(order)))
    stopifnot(abs(diff(match(c("A_near","Y_near"),order)))==1,
              abs(diff(match(c("B_far","Z_far"),order)))==1)
    stopifnot(!identical(order,sort(order)),identical(distances,original))
    built <- ggplot_build(plot)$data[[1]]
    # Matching sample identities lie on the top-left to bottom-right diagonal.
    zero <- built[built$fill=="#FFFFFF",,drop=FALSE]
    stopifnot(nrow(zero)==4,all(as.numeric(zero$x)+as.numeric(zero$y)==5))
    stopifnot(identical(plot$data$distance,as.vector(distances)))
    if (identical(permutation,1:4)) {
      stem <- file.path(output,paste0("clustered_",function_index))
      ggsave(paste0(stem,".png"),plot,width=8,height=7,dpi=120)
      ggsave(paste0(stem,".pdf"),plot,width=8,height=7)
    }
  }
}
cat("PASS: both methods cluster by VST distance; matching axes; top-left diagonal; input-order independence of group adjacency; unchanged numeric distances; PNG/PDF export.\n")
