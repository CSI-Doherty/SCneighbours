#' Check Single-Cell Object Type and Extract Components
#'
#' @title check_single_cell_object
#' @description Validates that the input object is either a Seurat object or a
#'   SingleCellExperiment object. Returns a named list containing standardized
#'   accessors for metadata, graphs, and other components needed by package functions.
#' @param obj An object to validate (should be either Seurat or SingleCellExperiment).
#' @return A named list with the following elements:
#'   \itemize{
#'     \item \code{type} - Character string: "Seurat" or "SingleCellExperiment"
#'     \item \code{embeddings} - Tibble of the selected reduction embeddings
#'     \item \code{key} - Character string of embedding dimension name key
#'     \item \code{metadata} - Data frame of cell metadata
#'     \item \code{graphs} - Seurat style nearest neighbour graph
#'     \item \code{n_cells} - Integer, number of cells in the dataset
#'     \item \code{cell_names} - Character vector of cell identifiers
#'   }
#' @export
#' @importFrom magrittr %>%
#' @importFrom igraph as_adjacency_matrix
#' @importFrom tibble as_tibble
#' @importFrom SingleCellExperiment reducedDim
check_single_cell_object <- function(obj, graph, reduction = NULL){
		if (inherits(obj, "scn_object")){
			return(obj)
		} else if (inherits(obj, "Seurat")) {
    	if(is.null(reduction)){reduction = 'umap'}
			if(inherits(graph, "character")){
				adjm = obj@graphs[[graph]]
				if(is.null(adjm)){
					stop(paste("Nearest Neighbour graph: ", graph, " not found.\n", "Please provide dgCMatrix, igraph, or select from graph names: ", paste(names(obj@graphs), collapse =  ", ")))
				}
			} else if(inherits(graph, "igraph")){
				adjm <- as_adjacency_matrix(graph)
				colnames(adjm) <- colnames(obj)
				rownames(adjm) <- colnames(obj)
			} else if(inherits(graph, "dgCMatrix")){
				adjm = graph
				colnames(adjm) <- colnames(obj)
				rownames(adjm) <- colnames(obj)
			} else {stop("Nearest neighbour either name of graph stored in Seurat object or adjacency graph class igraph or dgCMatrix must be included with SingleCellExperiment object. ",
										 "Received object of class: ", class(graph))}
			if(!(reduction %in% names(obj@reductions))){
				stop(paste("Reduction: ", reduction, "not found.\n"), paste("Please use one of the found reduction names: ", paste(names(obj@reductions), collapse = ', ')))
			}
        nc <-list(
            type = "Seurat",
            #object = obj,
            embeddings = obj@reductions[[reduction]]@cell.embeddings %>% as_tibble(rownames = "bc"),
            key = obj@reductions[[reduction]]@key,
            metadata = obj@meta.data,
            graph = adjm,
            n_cells = nrow(obj@meta.data),
            cell_names = colnames(obj)
        )
    } else if (inherits(obj, "SingleCellExperiment")) {
    	if(is.null(reduction)){reduction = 'UMAP'}
    	if(inherits(graph, "igraph")){
    		adjm <- as_adjacency_matrix(graph)
    		colnames(adjm) <- colnames(obj)
    		rownames(adjm) <- colnames(obj)
    	} else if(inherits(graph, "dgCMatrix")){
    		adjm = graph
    		colnames(adjm) <- colnames(obj)
    		rownames(adjm) <- colnames(obj)
    	} else {stop("Nearest neighbour adjacency graph of class igraph or dgCMatrix must be included with SingleCellExperiment object. ",
    							 "Received object of class: ", class(graph))}
    	nc <-list(
            type = "SingleCellExperiment",
            #object = obj,
            embeddings = reducedDim(obj, reduction) %>% as_tibble(rownames = "bc"),
            key = gsub('1', '',colnames(reducedDim(obj, reduction))[1]),
            metadata = as.data.frame(SingleCellExperiment::colData(obj)),
            graph = SeuratObject::as.Graph(adjm),
            n_cells = ncol(obj),
            cell_names = colnames(obj)
        )
    } else {
        stop("Input object must be either a Seurat object or a SingleCellExperiment object. ",
             "Received object of class: ", class(obj)[1])
    }
	class(nc) <- 'scn_object'
	return(nc)
}