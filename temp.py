### Create a heatmap of genes vs sample.
r("pdf('test.pdf')")
r("library(RColorBrewer)")
r("library(gplots)")
r("hmcol=colorRampPalette(brewer.pal(9,'GnBu'))(100)")
r("heatmap.2(exprs(vsd), col=hmcol, trace='none', margin=c(15,15), cexRow=0.8, cexCol=0.8)")