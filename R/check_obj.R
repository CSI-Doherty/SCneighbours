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
#'     \item \code{object} - The original input object
#'     \item \code{metadata} - Data frame of cell metadata
#'     \item \code{graphs} - List of available neighbor graphs
#'     \item \code{n_cells} - Integer, number of cells in the dataset
#'     \item \code{cell_names} - Character vector of cell identifiers
#'   }
#' @export
check_single_cell_object <- function(obj){

    if (inherits(obj, "Seurat")) {
        return(list(
            type = "Seurat",
            object = obj,
            metadata = obj@meta.data,
            graphs = obj@graphs,
            n_cells = nrow(obj@meta.data),
            cell_names = colnames(obj)
        ))
    } else if (inherits(obj, "SingleCellExperiment")) {
        return(list(
            type = "SingleCellExperiment",
            object = obj,
            metadata = as.data.frame(SingleCellExperiment::colData(obj)),
            graphs = SingleCellExperiment::colPairs(obj, asSparse = TRUE),
            n_cells = ncol(obj),
            cell_names = colnames(obj)
        ))
    } else {
        stop("Input object must be either a Seurat object or a SingleCellExperiment object. ",
             "Received object of class: ", class(obj)[1])
    }

}