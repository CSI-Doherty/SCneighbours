# visualise_neighbour_percentage

Creates a clustered heatmap showing the percentage of shared neighbours
between different cell types or groups. The heatmap is ordered by
hierarchical clustering to group similar neighbourhood patterns
together. This visualization helps identify which cell types have
similar spatial distributions or microenvironments.

## Usage

``` r
visualise_neighbour_percentage(obj, meta_data_column, graph = NULL)
```

## Arguments

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will first
  be converted to SCneighbours format

- meta_data_column:

  Name of the metadata column in the object to pull values from for
  identifying the cells of interest.

- graph:

  either a nearest neigbour graph in igraph, dgCMatrix or Seurat format,
  or the name of a graph stored in the Seurat object. (e.g., "RNA_nn",
  "RNA_snn", or "SCT_nn").

## Value

A ggplot2 object showing a heatmap of shared neighbour percentages
between cell types, with hierarchical clustering on both axes.

## Details

Visualize Neighbour Percentage Heatmap
