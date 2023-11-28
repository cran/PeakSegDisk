## ----setup, include = FALSE-----------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(width=100)

## -------------------------------------------------------------------------------------------------
data(Mono27ac, package="PeakSegDisk", envir=environment())
library(data.table)
loss.list <- list()
N.data.vec <- 10^seq(1, 3)
for(penalty in c(0, 1e2, 1e4, 1e6)){
  for(N.data in N.data.vec){
    some.cov <- Mono27ac$coverage[1:N.data]
    some.inc <- data.table(some.cov)
    some.inc[, count := 1:.N]
    data.list <- list(
      real=some.cov,
      synthetic=some.inc)
    for(short.type in names(data.list)){
      df <- data.list[[short.type]]
      fit <- PeakSegDisk::PeakSegFPOP_df(df, penalty)
      loss.list[[paste(penalty, N.data, short.type)]] <- data.table(
        N.data,
        short.type,
        fit$loss)
    }
  }
}
(loss <- do.call(rbind, loss.list))[, .(
  penalty, short.type, N.data, 
  mean.intervals, max.intervals,
  megabytes, seconds)][order(penalty, short.type, N.data)]

## -------------------------------------------------------------------------------------------------
(worst.dt <- data.table(
  N.data=N.data.vec,
  mean.intervals=(N.data.vec+1)/2,
  short.type="theoretical"))

## ----fig.height=3---------------------------------------------------------------------------------
one <- function(short.type, data.type, color){
  data.table(short.type, data.type, color)
}
type.dt <- rbind(
  one("theoretical", "Theoretical\nworst case", "grey"),
  one("synthetic", "Synthetic\nIncreasing", "red"),
  one("real", "Real ChIP-seq", "black"))
loss.types <- type.dt[loss, on=list(short.type)]
worst.types <- type.dt[worst.dt, on=list(short.type)]
(type.colors <- type.dt[, structure(color, names=data.type)])
if(require(ggplot2)){
ggplot()+
  guides(
    color=guide_legend(keyheight=3)
  )+
  geom_blank(aes(
    N.data, 1),
    data=data.table(N.data=c(5, 2000)))+
  theme_bw()+
  theme(panel.spacing=grid::unit(0, "lines"))+
  facet_grid(. ~ penalty, labeller=label_both)+
  geom_line(aes(
    N.data, mean.intervals, color=data.type),
    data=worst.types)+
  scale_color_manual(
    values=type.colors,
    breaks=names(type.colors))+
  geom_line(aes(
    bedGraph.lines, mean.intervals, color=data.type),
    data=loss.types)+
  scale_x_log10("N data")+
  scale_y_log10("Mean intervals (candidate changepoints)")
}

