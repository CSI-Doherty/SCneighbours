
#' @title neighbour_distance.scaled.
#' @description The function to calculate the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#' @param i The cell index.
#' @param reduction The reduction map used to calculate the coordinate variance.
#' @param seu The Seurat object.
#' @param graph Name of the nearest-neighbour graph to use from seu@graphs (default: "RNA_nn").
#' @return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#'
#' @export
#' @importFrom BBmisc normalize
#' @importFrom Seurat Embeddings
neighbour_distance.scaled = function(i, reduction, seu, graph = "RNA_nn") {
  # extract neighbour cells of cell i
  n = colnames(seu@graphs[[graph]])[seu@graphs[[graph]][i,] == 1]

  # normalize dimension scales
  embed.scale <- Seurat::Embeddings(seu, reduction = reduction)
  embed.scale[,1] <- BBmisc::normalize(embed.scale[,1], method = "range", range = c(-8, 8))
  embed.scale[,2] <- BBmisc::normalize(embed.scale[,2], method = "range", range = c(-8, 8))

  # dims of 20 cells
  e = embed.scale[n, ]

    if('numeric' %in% class(e)){
      
      return(mean(c(e[[1]], e[[2]])))
  }
  # variance in dim
  return(mean(var(e[,1]), var(e[,2])))
}


#' @title calculate_neighbour_distance_for_all_cells.
#' @description The function to calculate the average variance in coordinates of neighbour cells of all cells in the Seurat object based on the provided reduction map.
#' @param seu The Seurat object.
#' @param reduction The reduction map used to calculate the coordinate variance.
#' @param colname The column name to store the neighbourhood distance value in metadata.
#' @param graph Name of the nearest-neighbour graph to use from seu@graphs.
#' @return A Seurat object with a new metadata column storing the variance in coordinates of neighbour cells for each cell
#' @export
calculate_neighbour_distance_for_all_cells <- function(seu, reduction, colname, graph) {
  for(i in 1:nrow(seu@meta.data)){
      
      seu@meta.data[[colname]][i] = neighbour_distance.scaled(i, reduction, seu, graph)
      
      }
  seu
}
