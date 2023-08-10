library(gplots)
library(RColorBrewer)

#########################################################
### parse arguments
#########################################################
# args[1] <- counts matrix
# args[2] <- title
args <- commandArgs(trailingOnly = TRUE)

#########################################################
### reading in data and transform it to matrix format
#########################################################
data <- read.csv(args[1], sep="\t")
rnames <- data[,1]                            # assign labels in column 1 to "rnames"
mat_data <- data.matrix(data[,2:ncol(data)])  # transform column 2-ncol into a matrix
rownames(mat_data) <- rnames                  # assign row names

#########################################################
### customizing and plotting heatmap
#########################################################

pdf(file = paste(args[2], ".pdf", sep=""),
    width = 7,        # 5 x 300 pixels
    height = 7,
    #paper = "a4",
    pointsize = 8) # smaller font size

# changes the distance measure and clustering method
# NOTE: Matrix here not symmetrical. For symmetrical matrices
# only one distance and cluster could and SHOULD be defined.
# Distance options: euclidean (default), maximum, canberra, binary, minkowski, manhattan
# Cluster options: complete (default), single, average, mcquitty, median, centroid, ward
row_distance = dist(mat_data, method = "euclidean")
row_cluster = hclust(row_distance, method = "complete")
#col_distance = dist(t(mat_data), method = "manhattan")
#col_cluster = hclust(col_distance, method = "ward.D")

heatmap.2(mat_data,
  #cellnote = mat_data,  # same data set for cell labels
  main = args[2], # heat map title
  #notecol = "black",      # change font color of cell labels to black#
  density.info = "none",  # turns off density plot inside color legend
  trace = "none",         # turns off trace lines inside the heat map
  #col = bluered,
  col = "heat.colors",
  colsep = c(1:ncol(mat_data)), # line between each column
  rowsep = c(1:nrow(mat_data)), # line between each row
  sepwidth = c(0.01, 0.1), # line size between col, row
  margins=c(12,8), # margins for col, row labels
  offsetCol = 40, # push labels off plot by x char
  cexCol=1.4/log10(ncol(mat_data)),
  cexRow=1.4/log10(nrow(mat_data)),
  srtCol = -45 # rotate col labels
  )
  #breaks = seq(min(0.9, max(mat_data[!is.na(mat_data)]))-0.01,1,0.001))    # use on color palette defined earlier
  #Rowv = as.dendrogram(row_cluster), # apply default clustering method
  #Colv = as.dendrogram(col_cluster)) # apply default clustering method

#dev.off()
