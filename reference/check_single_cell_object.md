# check_single_cell_object

Takes a Seurat or SingleCellExperiment object and converts it into a
scn_object.

## Usage

``` r
check_single_cell_object(obj, graph = NULL, reduction = NULL)
```

## Arguments

- obj:

  An Seurat or SingleCellExperiment object to convert.

- graph:

  Name of the nearest-neighbour graph to use from seu@graphs (e.g.,
  "RNA_nn", "RNA_snn", or "SCT_nn").

- reduction:

  Name of the dimensionality reduction to use (default: "umap").

## Value

A S3 scn_object with the following elements:

- `type` - Character string: "Seurat" or "SingleCellExperiment"

- `embeddings` - Tibble of the selected reduction embeddings

- `key` - Character string of embedding dimension name key

- `metadata` - Data frame of cell metadata

- `graphs` - Seurat style nearest neighbour graph

- `n_cells` - Integer, number of cells in the dataset

- `cell_names` - Character vector of cell identifiers

## Details

Check Single-Cell Object Type and Extract Components
