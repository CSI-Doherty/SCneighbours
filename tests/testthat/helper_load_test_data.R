SeuratData::InstallData("ifnb")
SeuratData::InstallData("pbmc3k")


ifnb <- SeuratData::LoadData("ifnb")
pbmc <- SeuratData::LoadData("pbmc3k")

ifnb <- subset(ifnb, stim == 'CTRL')

seu <- merge(ifnb, pbmc)
seu <- Seurat::NormalizeData(seu)
seu <- Seurat::FindVariableFeatures(seu)
seu <- Seurat::ScaleData(seu)
seu <- Seurat::RunPCA(seu)

seu <- Seurat::FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- Seurat::FindClusters(seu, resolution = 0.5)


seu <- Seurat::RunUMAP(seu, dims = 1:30, reduction = "pca")
#Seurat::DimPlot(seu, group.by = c("orig.ident", "seurat_clusters"))

sce <- Seurat::as.SingleCellExperiment(SeuratObject::JoinLayers(seu))
sce.graph <- scran::buildSNNGraph(x = sce)

adj.matrix <- igraph::as_adjacency_matrix(sce.graph)

seu.graph <- seu@graphs$RNA_nn
