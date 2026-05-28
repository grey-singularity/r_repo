

#Download the data 
sdrf<-read.table("hypertrophy_sample_table.txt",header=TRUE,row.names = 1,as.is=TRUE,sep="\t") 

countdata<-read.table("hypertrophy_gene_counts.txt",header=TRUE,row.names = 1,as.is=TRUE,sep="\t") 

#Subset av prover
#countdata<-countdata[,16:29]
#sdrf<-sdrf[16:29,]

Treatment<-factor(sdrf$Treatment)
TimePoint<-factor(sdrf$TimePoint)


#filter edgeR
library(edgeR)
meanLog2CPM <- rowMeans(log2(cpm(countdata) + 1))
hist(meanLog2CPM)

sum(meanLog2CPM <= 1)

countdata <- countdata[meanLog2CPM > 1, ]
dim(countdata)

#Normalize deseq2

library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = countdata,
  design = ~ TimePoint+Treatment,
  colData = data.frame(
    TimePoint = TimePoint,Treatment=Treatment))

normCounts <- rlog(dds)
#Prepare count data for analysis
dge <- DGEList(countdata)
#Calculate normalization factors
dge <- calcNormFactors(dge)
#Define design matrix
designDF <- data.frame(TimePoint = TimePoint, Treatment=Treatment)

design <- model.matrix(~ 0 + TimePoint + Treatment, data = designDF) 
head(designDF)

head(design)
#Fit GLM
dge <- estimateDisp(dge, design, robust = TRUE)
fit <- glmQLFit(dge, design, robust = TRUE)
#Statistical testing
colnames(design)
contrastMat <- makeContrasts(TreatmentVsControl = TreatmentET,
 levels = design)
qlfRes <- glmQLFTest(fit, contrast = contrastMat)
topRes <- topTags(qlfRes, n = nrow(fit$counts))
topRes <- subset(topRes$table, abs(logFC) > 1 & FDR < 0.05)
print(head(topRes))
print(nrow(topRes))
# --- CONTINUING FROM YOUR PREVIOUS CODE ---

# 1. Extract the names of the significant differentially expressed genes
sig_genes <- rownames(topRes)
print(paste("Number of DEGs for clustering:", length(sig_genes)))

# 2. Extract the VST normalized counts for ONLY the significant genes
vst_matrix <- assay(normCounts)
sig_vst_matrix <- vst_matrix[sig_genes, ]
sig_vst_matrix <-scale(sig_vst_matrix)

set.seed(123)
# 3. Perform Hierarchical Clustering explicitly (Optional, to view the dendrogram of samples)
# We transpose the matrix t() to calculate the distance between SAMPLES, not genes.
tf.data<-t(sig_vst_matrix)
corD <-cor(tf.data, method = "spearman")
corD <- as.dist(1-corD)
sample_hclust <- hclust(corD, method = "average")
plot(sample_hclust)
sample_hclust<-as.dendrogram(sample_hclust)
plot(sample_hclust)
plot(sample_hclust, main = "Hierarchical Clustering of Samples (Spearman, Average)")
abline(h=0.5,col="red")
clusters <- cutree(as.hclust(sample_hclust),h=0.5)
print(clusters)
max(clusters)
#sample_dist <- dist(t(sig_vst_matrix), method = "euclidean")
#sample_hclust <- hclust(sample_dist, method = "complete")
#plot(sample_hclust, main = "Hierarchical Clustering of Samples (Euclidean, Complete)")

# to calculate the distance between genes., not SAMPLES
#sample_dist <- dist((sig_vst_matrix), method = "euclidean")
#sample_hclust <- hclust(sample_dist, method = "complete")
#plot(sample_hclust, main = "Hierarchical Clustering of Samples (Euclidean, Complete)")
# 4. Generate a Heatmap (This performs clustering automatically on both rows and columns)
# If you don't have pheatmap installed, run: install.packages("pheatmap")
library(pheatmap)

# Create an annotation dataframe to color-code the top of the heatmap
df_annotation <- data.frame(
  Treatment = Treatment,
  TimePoint = TimePoint
)
rownames(df_annotation) <- colnames(sig_vst_matrix)
breaks<-seq(-2,2,length.out=101)
# Draw the heatmap
pheatmap(sig_vst_matrix,
         cluster_rows = TRUE,       # Cluster the genes (rows)
         cluster_cols = TRUE,       # Cluster the samples (columns)
         annotation_col = df_annotation, 
         scale = "row",             # Scales gene expression to Z-scores to highlight relative changes
         show_rownames = FALSE, # Set to TRUE if you have <50 genes, otherwise it's unreadable
         breaks=breaks,
        main = "Heatmap of DEGs (ET vs Control)")
###############################################################################################
# --- ENHANCING THE GENE CLUSTERING ---

# 1. Run pheatmap, but this time save it to an object.
# We use 'cutree_rows = 2' to physically separate the genes into 2 main clusters 
# (typically representing the up-regulated and down-regulated modules).
hm_obj <- pheatmap(sig_vst_matrix,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   annotation_col = df_annotation,
                   scale = "row",
                   show_rownames = FALSE,
                   breaks=breaks,
                   cutree_rows = 9,  # Splits the gene dendrogram into 2 branches
                   silent = TRUE)    # Prevents it from plotting just yet

# 2. Extract the cluster assignments for each individual gene from the heatmap object
gene_clusters <- cutree(hm_obj$tree_row, k=9)

# 3. Create a row annotation dataframe so we can visualize these clusters
df_row_annotation <- data.frame(
  GeneModule = factor(paste("Cluster", gene_clusters))
)
# Ensure the rownames match the genes exactly
rownames(df_row_annotation) <- rownames(sig_vst_matrix) 

# 4. Redraw the final, enhanced heatmap with the new row annotations
pheatmap(sig_vst_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = df_annotation,
         annotation_row = df_row_annotation, # Adds the color bar for gene clusters on the left
         scale = "row",
         show_rownames = FALSE,
         breaks=breaks,
         cutree_rows = 9,
         main = "Enhanced Heatmap: ET vs Control with Gene Modules")

# 5. Extract the actual gene IDs for downstream biological analysis
cluster1_genes <- names(gene_clusters[gene_clusters == 1])
cluster2_genes <- names(gene_clusters[gene_clusters == 2])
cluster3_genes <- names(gene_clusters[gene_clusters == 3])
cluster4_genes <- names(gene_clusters[gene_clusters == 4])
cluster5_genes <- names(gene_clusters[gene_clusters == 5])
cluster6_genes <- names(gene_clusters[gene_clusters == 6])
cluster7_genes <- names(gene_clusters[gene_clusters == 7])
cluster8_genes <- names(gene_clusters[gene_clusters == 8])
cluster9_genes <- names(gene_clusters[gene_clusters == 9])

print(paste("Number of genes in Cluster 1:", length(cluster1_genes)))
print(paste("Number of genes in Cluster 2:", length(cluster2_genes)))
print(paste("Number of genes in Cluster 3:", length(cluster3_genes)))
print(paste("Number of genes in Cluster 4:", length(cluster4_genes)))
print(paste("Number of genes in Cluster 5:", length(cluster5_genes)))
print(paste("Number of genes in Cluster 6:", length(cluster6_genes)))
print(paste("Number of genes in Cluster 7:", length(cluster7_genes)))
print(paste("Number of genes in Cluster 8:", length(cluster8_genes)))
print(paste("Number of genes in Cluster 9:", length(cluster9_genes)))

# You can now view the Ensembl IDs for a specific module:
# head(cluster1_genes)
################################################################################################

################################################################################################
# --- ADDING K-MEANS CLUSTERING ---
# --- DETERMINING OPTIMAL k ---

# Install the required package if you haven't already by uncommenting the line below:
install.packages("rlang")
install.packages("factoextra")
library("rlang")
library("ggplot2")
library(factoextra)
packageVersion("rlang")
# 1. The Elbow Method (WCSS)
# This calculates the Within-Cluster Sum of Squares for k = 1 to 10
elbow_plot <- fviz_nbclust(sig_vst_matrix, kmeans, method = "wss") + 
  ggtitle("Elbow Method for Optimal k")

# Display the plot
print(elbow_plot)


# 2. The Silhouette Method
# This calculates the average silhouette width for k = 1 to 10
sil_plot <- fviz_nbclust(sig_vst_matrix, kmeans, method = "silhouette") + 
  ggtitle("Silhouette Score for Optimal k")

# Display the plot
print(sil_plot)
# 1. Set a random seed for reproducibility 
# (K-means starts with random centers, so a seed ensures you get the exact same result every time)
set.seed(123)

# 2. Define the number of clusters (k)
# We will use k = 4 here to look for slightly more granular regulatory patterns than the previous k=2.
k_clusters <- 5

# 3. Perform K-means clustering on the DEG matrix
km_res <- kmeans(sig_vst_matrix, centers = k_clusters, iter.max = 10000)

# 4. Extract the cluster assignments for each gene
kmeans_clusters <- km_res$cluster

# 5. Create a row annotation dataframe for the K-means clusters
df_kmeans_annotation <- data.frame(
  KMeans_Module = factor(paste("K-Cluster", kmeans_clusters))
)
rownames(df_kmeans_annotation) <- rownames(sig_vst_matrix)

# 6. Order the matrix so that genes in the same K-means cluster are visually grouped together
ordered_matrix <- sig_vst_matrix[order(kmeans_clusters), ]

# 7. Draw the heatmap mapped to K-means clusters
# NOTE: We set cluster_rows = FALSE because we want the rows ordered strictly by our K-means results, 
# not by a new hierarchical dendrogram.
pheatmap(ordered_matrix,
         cluster_rows = FALSE,                  # Turn off hierarchical clustering for genes
         cluster_cols = TRUE,                   # Keep hierarchical clustering for samples
         annotation_col = df_annotation,
         annotation_row = df_kmeans_annotation, # Show the K-means color bar on the left
         scale = "row",
         show_rownames = FALSE,
         breaks=breaks,
         main = paste("Heatmap of DEGs Grouped by K-means (k =", k_clusters, ")"))

# 8. Extract the actual gene IDs for each K-means cluster for downstream biological analysis
k_cluster1_genes <- names(kmeans_clusters[kmeans_clusters == 1])
k_cluster2_genes <- names(kmeans_clusters[kmeans_clusters == 2])
k_cluster3_genes <- names(kmeans_clusters[kmeans_clusters == 3])
k_cluster4_genes <- names(kmeans_clusters[kmeans_clusters == 4])
k_cluster5_genes <- names(kmeans_clusters[kmeans_clusters == 5])

print("--- K-Means Clustering Results ---")
print(paste("Number of genes in K-means Cluster 1:", length(k_cluster1_genes)))
print(paste("Number of genes in K-means Cluster 2:", length(k_cluster2_genes)))
print(paste("Number of genes in K-means Cluster 3:", length(k_cluster3_genes)))
print(paste("Number of genes in K-means Cluster 4:", length(k_cluster4_genes)))
print(paste("Number of genes in K-means Cluster 4:", length(k_cluster5_genes)))
# Example to view the Ensembl IDs for K-Cluster 1:
# head(k_cluster1_genes)


# 1. Count the number of genes in each K-means cluster
cluster_counts <- c(length(k_cluster1_genes), 
                    length(k_cluster2_genes), 
                    length(k_cluster3_genes), 
                    length(k_cluster4_genes), 
                    length(k_cluster5_genes))

# 2. Create text labels including the numbers of genes
pie_labels <- paste("Cluster", 1:5, "\n(", cluster_counts, " genes)")

# 3. Define colors for each cluster
cluster_colors <- c("#E41A1C",  # Cluster 1 = red
                    "#377EB8",  # Cluster 2 = blue
                    "#4DAF4A",  # Cluster 3 = green
                    "#984EA3",  # Cluster 4 = purple
                    "#FF7F00")  # Cluster 5 = orange

# 4. Draw the pie chart
pie(cluster_counts, 
    labels = pie_labels, 
    main = "Distribution of DEGs Across 5 K-means Clusters",
    col = cluster_colors)

# 5. Add a legend
legend("topright",
       legend = paste("Cluster", 1:5),
       fill = cluster_colors,
       title = "K-means Clusters",
       cex = 0.8,
       bty = "n")
pca<-prcomp(sig_vst_matrix)
pcadf<-data.frame(PC1=pca$x[,1],PC2=pca$x[,2],Cluster=as.factor(km_res$cluster))
#plot(pcadf,col=km$cluster,pch=19,xlab="PC1",ylab="PC2")
library(ggplot2)
PCAplot<-ggplot(pcadf,aes(PC1,PC2,color=Cluster))+geom_point(size=3)+theme_minimal()+ggtitle("K-means clustering")
PCAplot
################################################################################################

