#! ~/.conda/envs/Rdeseq2/bin/R

library("optparse")

# parses command-line arguments 
parse_opts <- function() {
    options = list(
        make_option("--matrix", dest="matrix", help="Deseq results with padj and log2FoldChange columns"),
        make_option("--alpha", type="numeric", default=0.05, dest="alpha", help="P-value and adjusted p-value significance threshold"),
        make_option("--lfc", type="numeric", default=1, dest="lfc", help="Log fold change threshold"),
        make_option("--valpha", type="numeric", default=0.4, dest="valpha", help="Point transparency"),
        make_option("--out", dest="out", help="output prefix for all files")
    )

    opts = parse_args(OptionParser(option_list=options))

    if (is.null(opts$matrix)) {
        stop("Input matrix file not specified.")
    }
    return(opts)
}

make_plot <- function(counts_matrix, lfc_thr=1, pv_thr=0.05, valpha=0.4, figname="volcano") {
    df <- read.csv(counts_matrix, sep=",", header=TRUE, row.names=1, stringsAsFactors=FALSE)

    # convert data into volcano plot arrays
    x <- as.numeric(df$log2FoldChange)
    y <- -log10(as.numeric(df$padj))

    # color logic
    colors <- rep("gray", nrow(df))
    #colors[y < -1*log10(pv_thr)] <- "gray"
    colors[x > lfc_thr & y > -1*log10(pv_thr)] <- "green"
    colors[x < -lfc_thr & y > -1*log10(pv_thr)] <- "red"

    # get x axis bounds
    xbound <- max(abs(range(x)))

    # scatter plot
    plot(x, y, col=colors, pch=16, cex=valpha, xlab="log2FoldChange", ylab="-log10(P)", xlim=c(-xbound, xbound))

    # draw volcano plot lines
    abline(h=-log10(pv_thr), col="gray", lty="dashed")
    abline(v=-lfc_thr, col="gray", lty="dashed")
    abline(v=lfc_thr, col="gray", lty="dashed")
}

main <- function(){
    options <- parse_opts()
    pdf(file = paste(options$out, ".pdf", sep=""), width = 7, height = 7, pointsize = 8)
    # Call the make_plot function with the parsed arguments
    make_plot(options$matrix, options$lfc, options$alpha, options$valpha, options$out)
}

main()

#dev.off()