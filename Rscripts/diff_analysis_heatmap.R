library(DESeq2)
library(ggplot2)
library(pheatmap)
library(RColorBrewer) # Assicurati che questa libreria sia caricata

counts_raw <- read.table("path/to/counts.txt", header=TRUE, row.names=1, comment.char="#")


count_cols <- grep("\\.bam$", colnames(counts_raw))
counts <- counts_raw[, count_cols]


colnames(counts) <- gsub("^X\\.", "", colnames(counts))  
colnames(counts) <- sub(".*UNIPD", "UNIPD", colnames(counts))
colnames(counts) <- gsub("\\.sorted\\.bam$", "", colnames(counts))
colnames(counts) <- gsub("\\.", "-", colnames(counts))


conditions <- read.table("path/to/conditions.txt", header=TRUE, row.names=1, stringsAsFactors=FALSE)
rownames(conditions) <- gsub("\\.", "-", rownames(conditions))

cat("Colnames(counts):\n"); print(colnames(counts))
cat("Rownames(conditions):\n"); print(rownames(conditions))
if (!all(colnames(counts) == rownames(conditions))) {
  stop("ERROR: Sample names in counts and conditions files do not match exactly!")
}

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = conditions,
                              design = ~ 1)

keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- estimateSizeFactors(dds)
vsd <- varianceStabilizingTransformation(dds, blind=TRUE)

#selection of top 3 variable genes
rv <- rowVars(assay(vsd))
top_genes <- head(order(rv, decreasing=TRUE), 30)

mat <- assay(vsd)[top_genes, ]

anno <- as.data.frame(colData(vsd)[, "condition", drop=FALSE])

pheatmap(mat,
         annotation_col=anno,
         show_rownames=TRUE,
         cluster_cols=TRUE,
         color = colorRampPalette(brewer.pal(n = 9, name = "Blues"))(100),
         filename = "heatmap_top30.png")
