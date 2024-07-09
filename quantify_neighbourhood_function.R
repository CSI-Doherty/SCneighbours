library(Seurat)
library(dplyr)
library(tidyverse)

# calculate the variance in coordinates of the neighbour cells for aech cell in a certain dimensionality reduction

#'@title neighbour_distance.scaled.
#'@description The function to calculate the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#'@param i The cell index.
#'@param reduction The reduction map used to calculate the coordinate variance.
#'@param seu The seurat object.
#'@return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
neighbour_distance.scaled = function(i, reduction, seu, graph = "RNA_nn") {
  # extract neighbour cells of cell i
  n = colnames(seu@graphs[[graph]])[seu@graphs[[graph]][i,] == 1] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform

  # normalize dimension scales
  embed.scale <- Embeddings(seu, reduction = reduction)
  embed.scale[,1] <- BBmisc::normalize(embed.scale[,1], method = "range", range = c(-8, 8))
  embed.scale[,2] <- BBmisc::normalize(embed.scale[,2], method = "range", range = c(-8, 8))

  # dims of 20 cells
  e = embed.scale[n, ]

  # variance in dim
  mean(var(e[,1]), var(e[,2]))
}



#'@title calculate_neighbour_distance_for_all_cells.
#'@description The function to calculate the average variance in coordinates of neighbour cells of all cells in the seurat object based on the provided reduction map.
#'@param colname The colname to store the neighbourhood distance value in metadata.
#'@param reduction The reduction map used to calculate the coordinate variance.
#'@param seu The seurat object.
#'@return A seurat object with a new metadata column storing the varaince in coodinates of neighbour cells for each cell
calculate_neighbour_distance_for_all_cells <- function(seu, reduction, colname, graph) {
  for(i in 1:nrow(seu@meta.data)) {seu@meta.data[[colname]][i] = neighbour_distance.scaled(i, reduction, seu, graph)}
  seu
}




