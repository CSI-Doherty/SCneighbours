# calculate_neighbour_distance_for_all_cells.

The function to calculate the average variance in coordinates of
neighbour cells of all cells in the Seurat object based on the provided
reduction map.

## Usage

``` r
calculate_neighbour_distance_for_all_cells(
  obj,
  reduction = NULL,
  colname,
  graph = NULL
)
```

## Arguments

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will
  firist be converted to SCneighbours format

- reduction:

  The reduction map used to calculate the coordinate variance.

- colname:

  The column name to store the neighbourhood distance value in metadata.

- graph:

  either a nearest neigbour graph in igraph, dgCMatrix or Seurat format,
  or the name of a graph stored in the Seurat object. (e.g., "RNA_nn",
  "RNA_snn", or "SCT_nn").

## Value

A Seurat or SingleCellExperiment object with a new metadata column
storing the variance in coordinates of neighbour cells for each cell
