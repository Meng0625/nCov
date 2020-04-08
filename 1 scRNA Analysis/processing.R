############### Step1. Installing R packages ###################
install.packages('Seurat')


############### Step2. Load in the datasets ###################
# Change your working directory
setwd('/path/to/work')
# Download the dataset from http://.
# Download the gene features from http://.
# Locate the dataset to the work path
count <- readRDS('/path/to/work/datasets.rds')
gene <- readRDS('/path/to/work/genes.rds')
rownames(count) <- gene[,2]
dim(count)
#23045 140956
count <- count[!duplicated(rownames(count)), ]
dim(count)
#23038 140956

############### Step3. Integration of multiple single-cell datasets ###################
library(Seurat)
library(dplyr)
library(pheatmap)
library(ggplot2)

# Initialize the Seurat object with the raw (non-normalized data)
seu <- CreateSeuratObject(
  count
)
seu@meta.data$sample<-sub('-.*','',colnames(count))
seu@meta.data$sample<-factor(seu@meta.data$sample,
                             levels = c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12',
                                        'HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8'),
                             ordered = TRUE)

# QC and selecting cells for further analysis
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

# Normalizing the data
seu <- NormalizeData(object = seu, normalization.method = "LogNormalize", scale.factor = 1e4)

# Identification of highly variable features (feature selection)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(seu), 10)

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(seu)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2


# To construct a reference, we will identify ‘anchors’ between the individual datasets. 
# First, we split the combined object into a list, with each dataset as an element.
seu.list <- SplitObject(seu, split.by = "sample")

seu.list <- seu.list[c('P5','P1','P2','P6','P3','P7','P8','P9','P10','P11','P4','P12',
                       'HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8')]

# Prior to finding anchors, we perform standard preprocessing (log-normalization), 
# and identify variable features individually for each.
for (i in 1:length(seu.list)) {
  seu.list[[i]] <- NormalizeData(seu.list[[i]], verbose = FALSE)
  seu.list[[i]] <- FindVariableFeatures(seu.list[[i]], selection.method = "vst", 
                                       nfeatures = 2000, verbose = FALSE)
}

# Identify anchors
seu.anchors <- FindIntegrationAnchors(object.list = seu.list, dims = 1:30)

# Integrated datasets
seu.integrated <- IntegrateData(anchorset = seu.anchors, dims = 1:30)


############### Step4. Cell Cycle Regression ##################
# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# Assign Cell-Cycle Scores
seu.integrated <- CellCycleScoring(seu.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Visualize the distribution of cell cycle markers across
RidgePlot(seu.integrated, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

# Regressing out the difference between the G2M and S phase scores.
seu.integrated$CC.Difference <- seu.integrated$S.Score - seu.integrated$G2M.Score
seu.integrated <- ScaleData(seu.integrated, vars.to.regress = "CC.Difference", features = rownames(seu), verbose = FALSE)

############### Step5. Clustering and visualization ##################
# Perform linear dimensional reduction
seu.integrated <- RunPCA(seu.integrated, npcs = 100, verbose = FALSE, ndims.print = 1:5, nfeatures.print = 5)
DimPlot(seu.integrated, reduction = "pca",label=TRUE)

# Generates an ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one
ElbowPlot(seu.integrated, ndims = 100)

# Heatmap for easy exploration of the primary sources of heterogeneity
DimHeatmap(seu.integrated, dims = c(15:20), cells = 500, balanced = TRUE)

# Clustering
# Construct a KNN graph based on the euclidean distance in PCA space, 
# and refine the edge weights between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity).
seu.integrated <- FindNeighbors(seu.integrated, dims = 1:30, nn.eps = 1)

# Cluster the cells
# resolution: Value of the resolution parameter, use a value above (below) 1.0 if you want to obtain a larger (smaller) number of communities.
seu.integrated <- FindClusters(seu.integrated, resolution = 1.2)

# Run non-linear dimensional reduction (tSNE)
seu.integrated <- RunTSNE(seu.integrated, dims = 1:20)
DimPlot(seu.integrated, reduction = "tsne",label=TRUE)

# Run non-linear dimensional reduction (UMAP)
seu.integrated <- RunUMAP(seu.integrated, reduction = "pca", dims = 1:20)
DimPlot(seu.integrated, reduction = "umap",label=TRUE)



############### Step6. Cluster biomarkers ##################
# Find markers for every cluster compared to all remaining cells
seu.integrated.markers <- FindAllMarkers(seu.integrated)

# Visualizing marker expression
VlnPlot(seu.integrated, features = c("CD4", "CD8A"))
FeaturePlot(seu.integrated, features = c("CD4", "CD8A"))



############### Step7. Assigning cell type identity to clusters ##################
seu.integrated$orig.ident<-seu.integrated$seurat_clusters

new.cluster.ids <- c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                     "14","15","16","17","18","19","20","21","22","23","24",
                     "25","N","26","27","N","N","28",
                     "N","N","N","N","N","N","N","N","N","N","N","N","N","N")

names(new.cluster.ids) <- levels(seu.integrated)
seu.integrated$new.cluster.ids<-new.cluster.ids[seu.integrated@active.ident]

new.cluster.ids.label <- c("C1: Naïve CD4+ T cells", 
                           "C2: NK cells", 
                           "C3: Monocytes", 
                           "C4: Monocytes", 
                           "C5: Monocytes", 
                           "C6: Monocytes", 
                           "C7: ADCC monocytes", 
                           "C8: CD4+ T cells",
                           "C9: NK cells", 
                           "C10: Naïve B cells", 
                           "C11: Cytotoxic CD8+ T cells", 
                           "C12: Monocytes",
                           "C13: CD4+ T cells", 
                           "C14: Monocytes",
                           "C15: Naïve CD8+ T cells",
                           "C16: Non-switched memory B cells",
                           "C17: Transitional CD8+ T cells",
                           "C18: Transitional CD8+ T cells",
                           "C19: Non classical monocytes (CD16+)",
                           "C20: CD4+ T cells",
                           "C21: Monocytes",
                           "C22: Intermediate monocytes",
                           "C23: Tregs",
                           "C24: Switched memory B cells",
                           "C25: CD8+ T cells",
                           "Others",
                           "C26: CD8+ T cells",
                           "C27: Plasma cells",
                           "Others",
                           "Others",
                           "C28: Dendritic cells",
                           "Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others")

names(new.cluster.ids.label) <- levels(seu.integrated)
seu.integrated$new.cluster.ids.label<-new.cluster.ids.label[seu.integrated@active.ident]

seu.integrated$new.cluster.ids<-factor(seu.integrated$new.cluster.ids,
                                       levels=c("1","2","3","4","5","6","7","8","9","10","11","12","13",
                                                "14","15","16","17","18","19","20","21","22","23","24",
                                                "25","26","27","28","N"),
                                       ordered = T)

bigger.clusters <- c("Conventional CD4+ T cells", #0 CD4+
                     "NK cells", #1
                     "Monocytes", #2
                     "Monocytes", #3
                     "Monocytes", #4
                     "Monocytes", #5
                     "Monocytes", #6
                     "Conventional CD4+ T cells", #7  CD4+
                     "NK cells", #8
                     "B cells", #9
                     "CD8+ T cells", #10 CD8+
                     "Monocytes", #11
                     "Conventional CD4+ T cells", #12  CD4+
                     "Monocytes", #13
                     "CD8+ T cells", #14 CD8+
                     "B cells", #15
                     "CD8+ T cells", #16 CD8+
                     "CD8+ T cells", #17 CD8+
                     "Monocytes", #18
                     "Conventional CD4+ T cells", #19  CD4+
                     "Monocytes", #20
                     "Monocytes", #21
                     "Regulatory CD4+ T cells", #22  Tregs
                     "B cells", #23
                     "CD8+ T cells", #24 CD8+
                     "Others", #25
                     "CD8+ T cells", #26 CD8+
                     "Plasma cells", #27
                     "Others", #28
                     "Others", #29
                     "Dendritic cells", #30
                     "Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others","Others")
names(bigger.clusters) <- levels(seu.integrated)
seu.integrated$bigger.cluster<-bigger.clusters[seu.integrated@active.ident]
seu.integrated$bigger.cluster<-factor(seu.integrated$bigger.cluster,
                                      levels=c("Monocytes",'Conventional CD4+ T cells','CD8+ T cells',
                                               'B cells','Plasma cells','NK cells','Dendritic cells', 'Regulatory CD4+ T cells'),
                                      ordered = T)

# Rename the clusters
seu.integrated <- RenameIdents(seu.integrated, new.cluster.ids.label)


# Add sample type
sampleType<-c('Normal','Normal','Normal','Normal','Normal','Normal','Normal','Normal',
              'COVID-19','COVID-19','COVID-19','COVID-19','COVID-19','COVID-19',
              'COVID-19','COVID-19','COVID-19','COVID-19','COVID-19','COVID-19')
names(sampleType)<-c('HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8',
                     'P1','P2','P3','P4','P5','P6',
                     'P7','P8','P9','P10','P11','P12')

seu.integrated$sampleType<-sampleType[seu.integrated$sample]
seu.integrated$sample<-factor(seu.integrated$sample,
                                   levels=c('P1','P2','P3','P4','P5','P6',
                                            'P7','P8','P9','P10','P11','P12',
                                            'HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8'),
                                   ordered=T)

# tSNE plot
DimPlot(seu.integrated, reduction = "tsne",label=TRUE)


############### Step8. Visualization ##################
cells<-colnames(seu.integrated)[seu.integrated$new.cluster.ids !='N']
cells_nCov<-colnames(seu.integrated)[seu.integrated$new.cluster.ids !='N' &
                                       seu.integrated$sampleType == 'COVID-19']
cells_norm<-colnames(seu.integrated)[seu.integrated$new.cluster.ids !='N' &
                                       seu.integrated$sampleType == 'Normal']

colors<-c('#5B9BD5', '#ED7D31', '#A5A5A5', '#FFC000', '#4472C4', '#70AD47', 
          '#255E91', '#9E480E', '#636363', '#99480E', '#264478', '#43682B', 
          '#7CAFDD', '#F1975A', '#B7B7B7', '#FFCD33', '#698ED0', '#8CC168', 
          '#327DC2', '#D26012', '#848484', '#CC9A00', '#335AA1', '#5A8A39', 
          '#9DC3E6', '#F4B183', '#C9C9C9', '#FFD966' )

##### Figure 1B
# pdf 15*8
DimPlot(seu.integrated, 
           reduction = "tsne",
           #group.by = 'new.cluster.ids',
           #cols=colors,
           #label = TRUE,
           cells = cells)

DimPlot(seu.integrated, reduction = "tsne", cells = cells_nCov)
DimPlot(seu.integrated, reduction = "tsne", cells = cells_norm)

##### Figure 1B (Cell Type)
# pdf 11x8
DimPlot(seu.integrated, 
        reduction = "tsne", 
        cells = cells, 
        group.by = 'bigger.cluster',
        #label.size = 6,
        #label=TRUE,
        cols=colors)
DimPlot(seu.integrated, 
        reduction = "tsne", 
        cells = cells_nCov, 
        group.by = 'bigger.cluster',
        label=TRUE,
        cols=colors)
DimPlot(seu.integrated, 
        reduction = "tsne", 
        cells = cells_norm, 
        group.by = 'bigger.cluster',
        label=TRUE,
        cols=colors)



##### Figure S1B
library(reshape2)
cluster_dis<-cbind(paste0('C',as.vector(seu.integrated$new.cluster.ids)),as.vector(seu.integrated$sample))
colnames(cluster_dis)<-c('cluster','sample')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=c('C3','C4','C5','C6','C7','C12','C14','C19','C21','C22',
                                     'C1','C8','C13','C20','C23','C11','C15','C17','C18',
                                     'C25','C26','C28','C2','C9','C10','C16','C24','C27','CN'),
                            ordered = TRUE)
cluster_dis$sample<-factor(cluster_dis$sample,
                           levels=c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12',
                                    'HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8'),
                           ordered = TRUE)
cluster_dis<-cluster_dis[cluster_dis$cluster != 'CN',]
ggplot(cluster_dis, aes(cluster)) + 
  geom_bar(aes(fill=sample), position="fill", color="black") +
  scale_fill_manual(values=c(rep('gray45',4),
                             rep('gray70',8),
                             rep('snow',8)))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


##### Figure 1E
library(reshape2)
cluster_dis<-cbind(as.vector(seu.integrated$sample),as.vector(seu.integrated$bigger.cluster))
colnames(cluster_dis)<-c('sample','cluster')
cluster_dis<-as.data.frame(cluster_dis)
cluster_dis<-cluster_dis[!is.na(cluster_dis$cluster),]
cluster_dis$cluster<-factor(cluster_dis$cluster,
                            levels=c('Monocytes','Conventional CD4+ T cells','CD8+ T cells','B cells',
                                     'Plasma cells','NK cells','Dendritic cells','Regulatory CD4+ T cells'),
                            ordered = TRUE)
cluster_dis$sample<-factor(cluster_dis$sample,
                           levels=c('P1','P2','P3','P4','P5','P6','P7','P8','P9','P10','P11','P12',
                                    'HC1','HC2','HC3','HC4','HC5','HC6','HC7','HC8'),
                           ordered = TRUE)

## pdf 8x3
ggplot(cluster_dis, aes(sample)) + 
  geom_bar(aes(fill=cluster), position="fill", color="black") +
  scale_fill_manual(values=colors)+
  labs(x = "Sample",
       y = "Percentage")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.title=element_blank())




##### Figure 1C
genes<-c("FOXP3","CTLA4","IL2RA",
             "CD7","CD2","CD3G","IL7R",
             "CD3E","CD3D","CCL5","NKG7","GNLY","PRF1","GZMA","CST7","KLRD1","KLRB1","SPON2",
             "KLRF1","NCR1","NCAM1","KIR2DL3",
             "CD163","VCAN","S100A12","CD14","SERPINA1","CD68","LST1","FCN1","CST3",
             #"LYZ","S100A9","S100A8",
             "FCGR2A","CLEC7A",
             "CLEC10A","CD1C","FCER1A","CD1E",
             "BANK1","MS4A1",
             "CD79A","IGKC","JCHAIN",
             "MZB1","IGHA1","IGLL5",
             "CD38","TNFRSF17")

seu.integrated_markers_m<-seu.integrated@assays$RNA@data[rownames(seu.integrated@assays$RNA@data) %in% genes,]
seu.integrated_markers_m<-t(as.matrix(seu.integrated_markers_m))
seu.integrated_markers_m<-as.data.frame(seu.integrated_markers_m)
seu.integrated_markers_m<-cbind(seu.integrated_markers_m,paste0('C',as.vector(seu.integrated$new.cluster.ids)))

colnames(seu.integrated_markers_m)[dim(seu.integrated_markers_m)[2]]<-'cluster'
#seu.integrated_markers_m$cluster<-as.vector(seu.integrated_markers_m$cluster)
seu.integrated_markers_m2<-aggregate(seu.integrated_markers_m[,1:(dim(seu.integrated_markers_m)[2]-1)],by=list(cluster=seu.integrated_markers_m$cluster),FUN=mean)
seu.integrated_markers_m3<-seu.integrated_markers_m2[,-1]

seu.integrated_markers_m4<-matrix(unlist(seu.integrated_markers_m3),dim(seu.integrated_markers_m3)[1],dim(seu.integrated_markers_m3)[2])
colnames(seu.integrated_markers_m4)<-colnames(seu.integrated_markers_m3)
rownames(seu.integrated_markers_m4)<-as.vector(seu.integrated_markers_m2[,1])
seu.integrated_markers_m4<-log2(seu.integrated_markers_m4+1)
pheatmap(t(seu.integrated_markers_m4),
         angle_col = 45)

seu.integrated_markers_m4_anno<-cbind(as.vector(seu.integrated$new.cluster.ids),
                                      as.vector(seu.integrated$new.cluster.ids.label),
                                      as.vector(seu.integrated$bigger.cluster))
rownames(seu.integrated_markers_m4_anno)<-paste0('C',seu.integrated_markers_m4_anno[,1])
colnames(seu.integrated_markers_m4_anno)<-c('clusterID','Cluster','Type')
seu.integrated_markers_m4_anno<-unique(seu.integrated_markers_m4_anno)
seu.integrated_markers_m4_anno<-as.data.frame(seu.integrated_markers_m4_anno)
seu.integrated_markers_m4_anno<-seu.integrated_markers_m4_anno[,-1]
seu.integrated_markers_m4_anno<-seu.integrated_markers_m4_anno[rownames(seu.integrated_markers_m4_anno)!='CN',]
seu.integrated_markers_m4_anno$Cluster<-factor(as.vector(seu.integrated_markers_m4_anno$Cluster))
pheat<-pheatmap(t(seu.integrated_markers_m4))
seu.integrated_markers_m4<-seu.integrated_markers_m4[,pheat$tree_row$order]
seu.integrated_markers_m4<-seu.integrated_markers_m4[pheat$tree_col$order,]
seu.integrated_markers_m4<-seu.integrated_markers_m4[rownames(seu.integrated_markers_m4)!='CN',]
#seu.integrated_markers_m4_anno$Cluster<-seu.integrated_markers_m4_anno$Type
anno_colors<-list(
  Type=c(colors[1:8]),
  #Cluster=c(rainbow(28,alpha=0.5))
  Cluster=c('#F8766D', '#EF7F46', '#E38900', '#D59100', '#C49A00', '#B0A100', '#99A800', 
            '#7CAE00', '#53B400', '#00B823', '#00BC56', '#00BF77', '#00C094', '#00C1AD', 
            '#00BFC4', '#00BCD9', '#00B6EB', '#00AEF9', '#06A4FF', '#7698FF', '#A58AFF', 
            '#C77CFF', '#DF70F8', '#F166E9', '#FB61D7', '#FF62C1', '#FF66A8', '#FE6D8C' )
)
names(anno_colors[['Type']])<-
  c('Monocytes','Conventional CD4+ T cells','CD8+ T cells','B cells',
    'Plasma cells','NK cells','Dendritic cells','Regulatory CD4+ T cells')
names(anno_colors[['Cluster']])<-
  c("C1: Naïve CD4+ T cells","C2: NK cells","C3: Monocytes","C4: Monocytes",
    "C5: Monocytes","C6: Monocytes","C7: ADCC monocytes","C8: CD4+ T cells","C9: NK cells",
    "C10: Naïve B cells","C11: Cytotoxic CD8+ T cells","C12: Monocytes", "C13: CD4+ T cells",
    "C14: Monocytes","C15: Naïve CD8+ T cells","C16: Non-switched memory B cells",
    "C17: Transitional CD8+ T cells","C18: Transitional CD8+ T cells",
    "C19: Non classical monocytes (CD16+)","C20: CD4+ T cells","C21: Monocytes", 
    "C22: Intermediate monocytes", "C23: Tregs","C24: Switched memory B cells","C25: CD8+ T cells",
    "C26: CD8+ T cells","C27: Plasma cells", "C28: Dendritic cells")
#unique(seu.integrated_markers_m4_anno$Cluster)

seu.integrated_markers_m4<-seu.integrated_markers_m4[,genes]

# change cluster order
cluster_list<-c("C20","C13","C1","C8","C23","C15","C18","C25","C17","C26","C11","C9",
                "C2","C22","C7","C5","C21",
                "C3","C12","C14","C4","C6","C19","C28","C10","C16","C24","C27")
seu.integrated_markers_m4<-seu.integrated_markers_m4[cluster_list,]

## pdf 10x12
pheatmap(t(seu.integrated_markers_m4),
         annotation_col = seu.integrated_markers_m4_anno,
         cluster_rows=F,
         cluster_cols=F,
         angle_col=45,
         annotation_colors = anno_colors,
         show_colnames =T
)


## Save file
saveRDS(seu.integrated,'seu.integrated.rds')




