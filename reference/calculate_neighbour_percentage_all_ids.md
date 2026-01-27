# calculate_neighbour_percentage_all_ids

Calculates the percentage composition of neighbours for all cell
types/groups in the dataset. This creates a matrix showing what
percentage of each cell type's neighbourhood is composed of each other
cell type, useful for understanding global neighbourhood composition
patterns.

## Usage

``` r
calculate_neighbour_percentage_all_ids(obj, meta_data_column, graph = "RNA_nn")
```

## Arguments

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will
  firist be converted to SCneighbours format

- meta_data_column:

  Name of the column in seu@meta.data to pull values from for grouping
  and analysis.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (e.g.,
  "RNA_nn", "RNA_snn", or "SCT_nn").

## Value

A data frame where the first column 'ids' contains all unique cell type
labels, and subsequent columns (named by cell type) contain the
percentage of neighbours belonging to each cell type.

## Details

Calculate Neighbour Percentages for All Cell Types
