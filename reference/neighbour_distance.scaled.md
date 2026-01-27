# neighbour_distance.scaled.

The function to calculate the average variance in coordinates of
neighbour cells of a certain cell i based on the provided reduction map.

## Usage

``` r
neighbour_distance.scaled(i, reduction, obj, graph = NULL)
```

## Arguments

- i:

  The cell index.

- reduction:

  The name of the reduction map stored in the object to be used to
  calculate the coordinate variance.

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will
  firist be converted to SCneighbours format

- graph:

  either a nearest neigbour graph in igraph, dgCMatrix or Seurat format,
  or the name of a graph stored in the Seurat object. (e.g., "RNA_nn",
  "RNA_snn", or "SCT_nn").

## Value

A number, the average variance in coordinates of neighbour cells of a
certain cell i based on the provided reduction map.
