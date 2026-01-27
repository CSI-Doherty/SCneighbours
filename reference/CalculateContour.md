# CalculateContour

Calculates kernel density estimate contour lines for a specific cell
population in a dimensionality reduction space. Uses 2D kernel density
estimation to compute contours at a specified percentile level, useful
for visualizing the spatial distribution of cell populations.

## Usage

``` r
CalculateContour(
  obj,
  meta_data_column,
  meta_data_highlight,
  reduction = "umap",
  percent = 95
)
```

## Arguments

- obj:

  A Seurat, SingleCellExperiment or SCNeighbours object containing
  single-cell data. If in Seurat or SingleCellExperiment form will first
  be converted to SCneighbours format

- meta_data_column:

  Name of the metadata column in the object to pull values from for
  identifying the cells of interest.

- meta_data_highlight:

  The specific value within meta_data_column to calculate contours for.

- reduction:

  Name of the dimensionality reduction to use (default: "umap").

- percent:

  Percentile level for the contour (default: 95). Higher values create
  tighter contours around the densest regions.

## Value

A data frame with x and y coordinates defining the contour line.

## Details

Calculate Density Contour Lines for Cell Populations
