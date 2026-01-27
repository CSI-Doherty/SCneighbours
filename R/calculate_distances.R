
#' @title neighbour_distance.scaled.
#' @description The function to calculate the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#' @param i The cell index.
#' @param reduction The name of the reduction map stored in the object to be used to calculate the coordinate variance.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param graph either a nearest neigbour graph in igraph, dgCMatrix or Seurat format, or the name of a graph stored in the Seurat object.
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#'
#' @export
#' @importFrom BBmisc normalize
#' @importFrom dplyr filter
#' @importFrom stats var
neighbour_distance.scaled = function(i, reduction, obj, graph = NULL) {
  # extract neighbour cells of cell i
	
	scn <- check_single_cell_object(obj, graph, reduction)
	
	g <- scn[['graph']]
	embed.scale <- scn[['embeddings']]
	
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
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param reduction The reduction map used to calculate the coordinate variance.
#' @param colname The column name to store the neighbourhood distance value in metadata.
#' @param graph either a nearest neigbour graph in igraph, dgCMatrix or Seurat format, or the name of a graph stored in the Seurat object.
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @return A Seurat or SingleCellExperiment object with a new metadata column storing the variance in coordinates of neighbour cells for each cell
#' @export
calculate_neighbour_distance_for_all_cells <- function(obj, reduction = NULL, colname, graph = NULL) {
	scn <- check_single_cell_object(obj, graph, reduction)
	
	meta <- scn[['metadata']]
	n <- scn[['n_cells']] 
	
  for(i in 1:n){
      
      meta[[colname]][i] = neighbour_distance.scaled(i, reduction, scn, graph)
      
  }
	if(inherits(obj, "Seurat")){
		obj@meta.data <- meta
		return(obj)
	} else if(inherits(obj, "SingleCellExperiment")){
		obj@colData <- S4Vectors::DataFrame(meta)
		return(obj)
	}
	return(meta)
}
