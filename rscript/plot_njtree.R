#!/usr/bin/env Rscript
library("ape")

# get the command line arguments
args <- commandArgs(trailingOnly = TRUE)
data_file=args[1]
treetype=args[2]
outgroup=args[3]
tree_file=args[4]
pdf_file=args[5]

# TODO remove when done testing
data_file <- "data/njtree/horse.data"
treetype <- "phylogram" # "phylogram" / "fan"
# outgroup <- "Esom_0226A"
outgroup <- "Przewalski_Paratype"
tree_file <- "data/njtree/horse-phylogram.tree"
pdf_file  <- "data/njtree/horse-phylogram.pdf"

# load the distance matrix
dist.mx <- as.matrix(read.table(data_file, head=T, check.names=FALSE))

# build the NJ tree
tr <- bionjs(dist.mx[,-c(1:2)])

# root the tree
tr <- root(tr, outgroup = outgroup, resolve.root = TRUE)

# resolve.root adds a zero-length branch below the MRCA of the ingroup, so lets
# find the outgroup branch and place the root in the middle of the branch
last <- length(tr$edge.length)
len <- tr$edge.length[last]

tr$edge.length[1] <- len/2     # the edge to the MRCA is always first
tr$edge.length[last] <- len/2  # the edge to the outgroup is always last

# sort the tree
tr <- ladderize(tr)

# save the tree data
write.tree(tr, file=tree_file)

# calcualte the plot size
numnodes <- length(dist.mx[,1])

if (treetype == "phylogram") {
    plot.height <- max(c(numnodes/5, 7))
    plot.width <- 10
} else {
    plot.height <- max(c(numnodes/10, 10))
    plot.width <- plot.height + 2
}

# get the sample metadata
meta <- read.table("sample-metadata.tsv", sep="\t", quote='', header=TRUE, comment.char="")

# join the metadata, so we can get the node colours
info <- merge(dist.mx, meta[c('Sample','Type.Name','Colour','Order')], by = "Sample", all.x=TRUE)
rownames(info) <- info[,'Sample']

# get the population names and colours for the legend
key <- unique(info[c('Type.Name','Colour','Order')])
key <- key[with(key, order(Order)), ]
key[] <- lapply(key, as.character)

# fix the type issue with the colour codes
cols <- sapply(info[tr$tip.label, 'Colour'], as.character)

pdf(file=pdf_file, width = plot.width, height = plot.height)

# plot the tree
plot(tr, type=treetype, tip.color=cols, label.offset=0.001)
legend(x="topright",legend = key[[1]], fill = key[[2]], cex=0.8)
add.scale.bar()

dev.off()
