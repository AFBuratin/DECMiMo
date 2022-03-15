library(RColorBrewer)
# display.brewer.all()

cols <- c(
  # DEseq
  brewer.pal(n = 9, "YlOrRd")[c(5)],
  # Edger
  brewer.pal(n = 9, "GnBu")[c(7)],
  # limma
  brewer.pal(n = 9, "RdPu")[c(5)]
)

methods <- c("DESeq2",
             "edgeR",
             "voom")

names(cols) <- methods

methods2 <- c("DESeq2",
              "edgeR",
              "voom")

names(cols) <- methods2

method_cols <- c(
  brewer.pal(n=4, "Dark2")
)

detection_methods <- c("circexplorer2",
                       "ciri",
                       "dcc",
                       "findcirc"
                       )
names(method_cols) <- detection_methods

