# calculate_outside_neighbours_cell

For each cell, calculates what percentage of its neighbours belong to a
different cell type/group than itself. This metric quantifies cellular
heterogeneity in the local neighbourhood and can identify boundary cells
or cells in mixed microenvironments. The result is stored as a new
column in the Seurat object's metadata.

## Usage

``` r
calculate_outside_neighbours_cell(
  obj,
  meta_data_column,
  graph = "RNA_nn",
  colname
)
```

## Arguments

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will
  firist be converted to SCneighbours format

- meta_data_column:

  Name of the metadata column in the object to pull values from for
  grouping and analysis.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (e.g.,
  "RNA_nn", "RNA_snn", or "SCT_nn").

- colname:

  Name of the new metadata column to store the calculated percentage of
  outside neighbours for each cell.

## Value

The Seurat object with a new column in seu@meta.data containing the
percentage of neighbours (0-100) that belong to a different group than
the cell itself.

## Details

Calculate Percentage of Neighbours Outside Cell's Own Group
