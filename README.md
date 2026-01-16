# SCneighbours
Single-Cell Neighbourhood Exploration Tool 

Note: The functions are using the nearest-neighbour graph from the RNA assay (RNA_nn), PLEASE change to SCT_nn if you are using the SCT assay.

## Installation

library(devtools)

install_github("https://github.com/CSI-Doherty/SCneighbours")

## Technique 1: Visualization of Cell Neighbourhoods

For individual cells or cell groups, we visualize their cell neighbourhoods as either dots or density contours on a dimensionality reduction map, uncovering tightly interwoven or widely dispersed regions. 

This allows us to learn which clusters might need merging if they have largely overlapped neighbourhoods, and infer trajectory relationships of cells obscured by the reduction.
```
source("visualize_neighbourhood_function.R")

# example usage
visualize_neighbourhood(seu, meta_data_column = "seurat_clusters", meta_data_highlight = "12", reduction = "umap", density = T)
```

## Technique 2: Quantification of Cell Neighbourhoods

For each cell, we quantify the distribution of its neighbouring cells by the variance in their coordinates, providing insights into the cell positioning in the reduction map. 

Large variance indicates diminished confidence in a cell’s positioning.

```
source("quantify_neighbourhood_function.R")

# example usage
calculate_neighbour_distance_for_all_cells(seu, reduction = "umap", colname = "neighbourhood_variance")
```

## Technique 3: Complementary reductions

The cell neighbourhoods can be explored across various dimensionality reductions for insights. 

For example, diffusion map and force-directed layout can provide clues about cell development by preserving the continuous cell transitioning trajectories.


## Example

Here we use 2 different PBMC datasets from SeuratData to demonstrate usage.

```
SeuratData::InstallData("ifnb")
SeuratData::InstallData("pbmc3k")

ifnb <- SeuratData::LoadData("ifnb")
pbmc <- SeuratData::LoadData("pbmc3k")

ifnb <- subset(ifnb, stim == 'CTRL')
```

The PBMCs are merged without integration and Run through a standard pipeline to generate a nearest neighbour graph, clusters, and UMAP reduction.

```
seu <- merge(ifnb, pbmc)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu)
seu <- ScaleData(seu)
seu <- RunPCA(seu)

seu <- FindNeighbors(seu, dims = 1:30, reduction = "pca")
seu <- FindClusters(seu, resolution = 0.5)


seu <- RunUMAP(seu, dims = 1:30, reduction = "pca")
```

![Example UMAP](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/5a50c439-2851-4566-b4f6-8d889b42e404.png)
![CD3 UMAP](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/6d470ce6-5f3b-45d7-aac5-1104a38306c3.png)

```
DimPlot(seu, group.by = c("orig.ident", "seurat_clusters"))
FeaturePlot(seu, "CD3E")
```

Here the T cells are split into 2 clusters for the seperate datasets.

Using SC neighbours we can see that there is neighbourhood sharing between the 2 seperate T cell clusters. 
```
visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 1, 'umap', density = T) +
visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 1, 'umap', density = F)

visualize_neighbourhood(seu, 'seurat_clusters', 2, density = T, percent = 90) +
visualize_neighbourhood(seu, 'seurat_clusters', 2, density = F)
```
![Cluster 1 neighbourhoood](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/f3d27916-5923-4a63-81f6-4f2ca9c5587b.png)
![Cluster 2 neighbourhoood](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/d4e4f029-6729-40d5-abba-015398fb2b0a.png)

A heatmap showing all neighbourhood sharing of clusters can be generated with visualise_neighbour_percentage.
```
visualise_neighbour_percentage(seu = seu, graph = 'RNA_nn', meta_data_column = 'seurat_clusters')
```
![Heatmap Example](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/90951f5b-ea10-4cf6-8c51-60d88cd4be62.png)

For each cell we can calculate the percentage of cells from its neighbours share the cells cluster (or other metadata column).
```
seu <- calculate_outside_neighbours_cell(seu, 'seurat_clusters',colname = 'outside_neighbourhood')
FeaturePlot(seu, outside_neighbourhood)
```
![Outside neighbours](https://github.com/CSI-Doherty/SCneighbours/blob/cc94e9650441d4a5b53e13fd94244557d5b750de/img/823f2a7f-3fbc-477d-8be8-e9a472dae282.png)

Here we see can see the high percentage cells form the boundaries between clusters, even when those boundaries don't appear close together on the UMAP.

For each cell we can calculate the variance in the distance to it it neighbours based on a reductions embeddings.  This can be useful to see if there are any problems with your UMAP or there are some problamtic cells in your dataset.

```
seu <- calculate_neighbour_distance_for_all_cells(seu, reduction = 'umap', colname = 'neighbour_distance')
FeaturePlot(seu, 'neighbour_distance')
```



