### Named list of character vectors (column names of bed/bedGraph/tsv
### files), used to read data files, which do not contain a header /
### column names. Each name corresponds to a data/file type, and each
### value is a character vector of column names expected in that
### file. loss is for the coverage.bedGraph_penalty=VALUE_loss.tsv
### file generated by PeakSegFPOP_file; segments is for the
### coverage.bedGraph_penalty=VALUE_segments.bed generated by
### PeakSegFPOP_file; coverage is for the coverage.bedGraph file which
### is used as input to PeakSegFPOP_file.
col.name.list <- structure(
  list(
    loss=c(
      "penalty", "segments", "peaks", "bases", "bedGraph.lines",
      "mean.pen.cost", "total.loss", "equality.constraints",
      "mean.intervals", "max.intervals"),
    segments=c("chrom","chromStart", "chromEnd", "status", "mean"),
    coverage=c("chrom", "chromStart", "chromEnd", "count")
  ), ex=function(){

    library(PeakSegDisk)
    r <- function(chrom, chromStart, chromEnd, coverage){
      data.frame(chrom, chromStart, chromEnd, coverage)
    }
    four <- rbind(
      r("chr1", 0, 10,  2),
      r("chr1", 10, 20, 10),
      r("chr1", 20, 30, 14),
      r("chr1", 30, 40, 13))
    write.table(
      four, tmp <- tempfile(),
      sep="\t", row.names=FALSE, col.names=FALSE)
    read.table(
      tmp, col.names=col.name.list$coverage)
    
    pstr <- "10.5"
    PeakSegFPOP_file(tmp, pstr)
    outf <- function(suffix){
      paste0(tmp, "_penalty=", pstr, "_", suffix)
    }
    fread.first(outf("segments.bed"), col.name.list$segments)
    fread.first(outf("loss.tsv"), col.name.list$loss)
    
  })

