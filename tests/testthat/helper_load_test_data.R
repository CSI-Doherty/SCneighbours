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
seu <- Seurat::FindClusters(seu, resolution = 0.5, cluster.name = "unintegrated_clusters")


seu <- Seurat::RunUMAP(seu, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
#Seurat::DimPlot(seu, reduction = "umap.unintegrated", group.by = c("stim", "unintegrated_clusters"), label = T)

sce <- Seurat::as.SingleCellExperiment(SeuratObject::JoinLayers(seu))
sce.graph <- scran::buildSNNGraph(x = sce, )