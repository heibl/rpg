tr <- rtree(12)

## visualize
plot(tr, type = "unrooted", no.margin = TRUE)
edgelabels(cex = 0.5, frame = "n", col = "red")
nodelabels(cex = 0.5, frame = "n", col = "blue")
