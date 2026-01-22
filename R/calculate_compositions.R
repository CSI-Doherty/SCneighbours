#' Calculate Neighbour Cell Type Percentages
#'
#' @title calculate_neighbour_percentage
#' @description Calculates the percentage composition of cell types among the
#'   neighbours of cells belonging to a specific group. Returns a data frame
#'   showing what percentage of the neighbourhood belongs to each cell type.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param meta_data_column Name of the column in seu@meta.data to pull values from
#'   for grouping and analysis.
#' @param meta_data_highlight The specific value within meta_data_column to
#'   highlight and analyze neighbours for.
#' @param graph either a nearest neigbour graph in igraph, dgCMatrix or Seurat format, or the name of a graph stored in the Seurat object.
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @return A data frame with columns 'ids' (cell type labels), 'Freq' (count),
#'   and 'f' (percentage of total neighbours).
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate
#' @export
calculate_neighbour_percentage <- function(obj, meta_data_column, meta_data_highlight, graph = 'RNA_nn'){
    scn = check_single_cell_object(obj, graph, reduction = NULL)
    
    g = scn[['graph']]
    meta = scn[['metadata']]
    
    sn = colnames(g)[Matrix::colSums(g[meta[[meta_data_column]] == meta_data_highlight,]) > 0]
    ids = factor(meta[sn,meta_data_column], levels = levels(factor(meta[,meta_data_column])))
    table(ids) %>% as.data.frame() %>%
        mutate(f = Freq/sum(Freq)*100)
}

#' Calculate Neighbour Percentages for All Cell Types
#'
#' @title calculate_neighbour_percentage_all_ids
#' @description Calculates the percentage composition of neighbours for all cell
#'   types/groups in the dataset. This creates a matrix showing what percentage
#'   of each cell type's neighbourhood is composed of each other cell type,
#'   useful for understanding global neighbourhood composition patterns.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param meta_data_column Name of the column in seu@meta.data to pull values from
#'   for grouping and analysis.
#' @param graph Name of the nearest-neighbour graph to use from seu@graphs
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @return A data frame where the first column 'ids' contains all unique cell type
#'   labels, and subsequent columns (named by cell type) contain the percentage
#'   of neighbours belonging to each cell type.
#' @export
calculate_neighbour_percentage_all_ids <- function(obj, meta_data_column, graph = 'RNA_nn'){
	scn = check_single_cell_object(obj, graph, reduction = NULL)
	
	g = scn[['graph']]
	meta = scn[['metadata']]
	
	
	ids = levels(factor(meta[,meta_data_column]))
	results = data.frame(ids = ids)
	for(i in ids){
		results[,i] = calculate_neighbour_percentage(scn, meta_data_column = meta_data_column, meta_data_highlight = i, graph)$f
	}
	return(results)
}



#' Calculate Percentage of Neighbours Outside Cell's Own Group
#'
#' @title calculate_outside_neighbours_cell
#' @description For each cell, calculates what percentage of its neighbours belong
#'   to a different cell type/group than itself. This metric quantifies cellular
#'   heterogeneity in the local neighbourhood and can identify boundary cells or
#'   cells in mixed microenvironments. The result is stored as a new column in
#'   the Seurat object's metadata.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param meta_data_column Name of the metadata column in the object to pull values from
#'   for grouping and analysis.
#' @param graph Name of the nearest-neighbour graph to use from seu@graphs
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @param colname Name of the new metadata column to store the calculated
#'   percentage of outside neighbours for each cell.
#' @return The Seurat object with a new column in seu@meta.data containing the
#'   percentage of neighbours (0-100) that belong to a different group than the
#'   cell itself.
#' @export
calculate_outside_neighbours_cell <- function(obj, meta_data_column, graph = "RNA_nn", colname){
	scn = check_single_cell_object(obj, graph, reduction = NULL)
	
	g = scn[['graph']]
	meta = scn[['metadata']]  
	
	
	# for(i in 1:nrow(scn@meta.data)){
	# 	ids = scn@meta.data[[meta_data_column]][scn@graphs[[graph]][i,]>0]
	# 	scn@meta.data[[colname]][i] = (1-sum(ids == scn@meta.data[[meta_data_column]][i])/length(ids))*100
	# }
	
	for(i in 1:nrow(meta)){
		ids = meta[[meta_data_column]][g[i,]>0]
		if(inherits(obj, 'Seurat')){
			obj@meta.data[[colname]][i] = (1-sum(ids == meta[[meta_data_column]][i])/length(ids))*100
		} else {
			obj[[colname]][i] = (1-sum(ids == meta[[meta_data_column]][i])/length(ids))*100
		}
		
	}
	return(obj)
}

