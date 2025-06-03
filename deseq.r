

library("DESeq2")
library(data.table)
library(apeglm)
library(ggplot2)
library("RColorBrewer")
library("pheatmap")

# parse arguments:
# args[1] <- counts csv
# args[2] <- affected prefix
# args[3] <- alpha value
args <- commandArgs(trailingOnly = TRUE)

if ( length(args) < 3 ) {
    print("Provide counts csv, affected prefix, and alpha value")
    quit()
}

# read counts into matrix from csv
counts <- as.matrix(read.csv(args[1], sep="\t", row.names=1))
counts <- round(counts, 0)
colnames(counts) <- sub("X", "", colnames(counts))

# extract sample names to create coldata matrix
samples <- c(colnames(counts))
conditions <- c()

# loop for separating affected and unaffected samples
i = 1
for (sample in samples) {
    condition <- "unaffected"
    if (grepl(args[2], sample, fixed=FALSE)) {
        condition <- "affected"
    }

    conditions[[i]] <- condition
    i <- i + 1
}

coldata <- data.frame(as.matrix(data.table(
    #"samples" = samples,
    "condition" = conditions
)))

coldata <- as.data.frame(lapply(coldata, unlist))
coldata$condition <- factor(coldata$condition)
rownames(coldata) <- samples

# print stats and checks
#print(head(counts,2))
#print(coldata)
#print(ncol(counts))
#print(nrow(coldata))
#print(class(counts))
#print(class(coldata))
#print(all(rownames(coldata) %in% colnames(counts)))

# performing deseq
dds <- DESeqDataSetFromMatrix(counts, DataFrame(coldata), ~ condition)
dds$condition <- relevel(dds$condition, ref = "unaffected")
dds <- DESeq(dds)

# generate stats summary and perform lfc shrinkage
res05 <- results(dds, alpha=args[3])
resLFC <- lfcShrink(dds, coef=2, res=res05, type="apeglm")

# order results by smallest padj
res05 <- res05[order(res05$padj),]
resLFC <- resLFC[order(resLFC$padj),]

# MA plot
plotMA(res05, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))

# isolate up and downregulated genes using adjusted pvalue
get_significant_genes <- function(df, pvalue = 0.05, l2fc = 1) {
  sig <- df[df$padj <= pvalue, ]
  up <- sig[sig$log2FoldChange >= l2fc, ]
  down <- sig[sig$log2FoldChange <= -l2fc, ]
  return(list("up"=up, "down"=down))
}

# writing out csvs
write.csv(res05, file="raw_results.csv", quote=FALSE)
write.csv(na.omit(res05), file="deseq_results_non_nan.csv", quote=FALSE)
write.csv(na.omit(resLFC), file="shrunk_deseq_results_non_nan.csv", quote=FALSE)

sig_genes <- get_significant_genes(na.omit(res05), args[3], 1)
shrunk_sig_genes <- get_significant_genes(na.omit(resLFC), args[3], 1)

write.csv(sig_genes$up, file="upregulated.csv", quote=FALSE)
write.csv(sig_genes$down, file="downregulated.csv", quote=FALSE)
write.csv(shrunk_sig_genes$up, file="shrunk_upregulated.csv", quote=FALSE)
write.csv(shrunk_sig_genes$down, file="shrunk_downregulated.csv", quote=FALSE)

# variance stabilizing transformation
vsd <- vst(dds, blind=TRUE)

# Create heat map from data
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- samples
colnames(sampleDistMatrix) <- samples
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
heatMap <- pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}

save_pheatmap_pdf(heatMap, "heatmap.pdf")

# PCA graphing function
pcaData <- plotPCA(vsd, returnData=TRUE)

# get percent variation
percentVar <- round(100 * attr(pcaData, "percentVar"))

# pca code
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = condition, name = name)) +
    stat_ellipse(aes(group = condition, linetype = condition), type = "t", level = 0.95, size = 1.25, show.legend = FALSE) +
    geom_point(size = 3, show.legend = TRUE) + 
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    coord_fixed() + theme_bw() + ggtitle("PCA") +
    guides(shape = guide_legend(order = 1),color = guide_legend(order = 2)) 

ggsave("PCA.png", width = 8, height = 4.6, dpi = 600)
