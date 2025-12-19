
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
#' @importFrom dplyr filter
neighbour_distance.scaled = function(i, reduction, seu, graph = NULL) {
  # extract neighbour cells of cell i
	
	obj <- check_single_cell_object(seu, graph, reduction)
	
	g <- obj[['graph']]
	embed.scale <- obj[['embeddings']]
	
  n = colnames(g)[g[i,] == 1]

  # normalize dimension scales
  embed.scale[,2] <- BBmisc::normalize(embed.scale[,2], method = "range", range = c(-8, 8))
  embed.scale[,3] <- BBmisc::normalize(embed.scale[,3], method = "range", range = c(-8, 8))

  # dims of 20 cells
  e <- embed.scale %>% filter(bc %in% n)

  if('numeric' %in% class(e)){
  	return(mean(c(e[[2]], e[[3]])))
  }
  # variance in dim
  return(mean(var(e[,2]), var(e[,3])))
}


#' @title calculate_neighbour_distance_for_all_cells.
#' @description The function to calculate the average variance in coordinates of neighbour cells of all cells in the Seurat object based on the provided reduction map.
#' @param seu The Seurat object.
#' @param reduction The reduction map used to calculate the coordinate variance.
#' @param colname The column name to store the neighbourhood distance value in metadata.
#' @param graph Name of the nearest-neighbour graph to use from seu@graphs.
#' @return A Seurat or SingleCellExperiment object with a new metadata column storing the variance in coordinates of neighbour cells for each cell
#' @export
calculate_neighbour_distance_for_all_cells <- function(seu, reduction, colname, graph = NULL) {
	obj <- check_single_cell_object(seu, graph, reduction)
	
	meta <- obj[['metadata']]
	n <- obj[['n_cells']] 
	
  for(i in 1:n){
      
      meta[[colname]][i] = neighbour_distance.scaled(i, reduction, obj, graph)
      
  }
	if(inherits(seu, "Seurat")){
		seu@meta.data <- meta
		return(seu)
	} else if(inherits(seu, "SingleCellExperiment")){
		seu@colData <- S4Vectors::DataFrame(meta)
		return(seu)
	}
	return(meta)
}
