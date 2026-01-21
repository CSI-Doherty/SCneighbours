#' Visualize Neighbour Percentage Heatmap
#'
#' @title visualise_neighbour_percentage
#' @description Creates a clustered heatmap showing the percentage of shared
#'   neighbours between different cell types or groups. The heatmap is ordered
#'   by hierarchical clustering to group similar neighbourhood patterns together.
#'   This visualization helps identify which cell types have similar spatial
#'   distributions or microenvironments.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will first be converted to SCneighbours format
#' @param meta_data_column Name of the metadata column in the object to pull values from
#'   for identifying the cells of interest.
#' @param graph either a nearest neigbour graph in igraph, dgCMatrix or Seurat format, or the name of a graph stored in the Seurat object.
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @return A ggplot2 object showing a heatmap of shared neighbour percentages
#'   between cell types, with hierarchical clustering on both axes.
#' @export
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_tile scale_fill_gradientn theme_classic theme element_text labs element_blank
#' @importFrom stats dist hclust
#' @importFrom magrittr %>%
visualise_neighbour_percentage <- function(obj, meta_data_column, graph = NULL) {
  scn <- check_single_cell_object(scn, graph, reduction)
  
	x = calculate_neighbour_percentage_all_ids(scn, meta_data_column, graph)
	d = dist(t(x[,-1]))
	h = hclust(d)

	x %>% pivot_longer(-ids) %>%
		mutate(name = factor(name, levels = h$labels[h$order]),
					 ids = factor(ids, levels = h$labels[h$order])) %>%
		ggplot(aes(name,ids)) +
		geom_tile(aes(fill = value)) +
		scale_fill_gradientn(colours = c("blue", "white", "red")) +
		theme_classic() +
		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
		#labs(x = element_blank(), y=element_blank(), fill = "% shared")
		labs(x = NULL, y= NULL, fill = "% shared")
}

#' Calculate Density Contour Lines for Cell Populations
#'
#' @title CalculateContour
#' @description Calculates kernel density estimate contour lines for a specific
#'   cell population in a dimensionality reduction space. Uses 2D kernel density
#'   estimation to compute contours at a specified percentile level, useful for
#'   visualizing the spatial distribution of cell populations.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will first be converted to SCneighbours format
#' @param meta_data_column Name of the metadata column in the object to pull values from
#'   for identifying the cells of interest.
#' @param meta_data_highlight The specific value within meta_data_column to
#'   calculate contours for.
#' @param reduction Name of the dimensionality reduction to use (default: "umap").
#' @param percent Percentile level for the contour (default: 95). Higher values
#'   create tighter contours around the densest regions.
#' @return A data frame with x and y coordinates defining the contour line.
#' @importFrom Seurat Embeddings
#' @importFrom ks kde
#' @importFrom grDevices contourLines
#' @importFrom magrittr %>%
CalculateContour = function(scn, meta_data_column, meta_data_highlight, reduction = "umap", percent = 95) {
	# d = Embeddings(scn, reduction = reduction) %>% as_tibble(rownames = "bc")

	obj <- check_single_cell_object(scn, graph, reduction)
	
	d <- obj[['embeddings']] 
	meta <- obj[['metadata']]
	
  #kd <- ks::kde(d[scn[[meta_data_column]] == meta_data_highlight,2:3], compute.cont=TRUE)
	
	x <- d[meta[[meta_data_column]] == meta_data_highlight,2:3]
	if(nrow(x) < 3){
		stop("Number of highlighted cell must be > 2")
	}
	
	kd <- ks::kde(x, compute.cont=TRUE)
	contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                      z=estimate, levels=cont[paste0(100-percent,"%")])[[1]])
  contour_95 <- data.frame(contour_95)
  #
  # d[scn[[meta_data_column]] == meta_data_highlight,2:3]
  return(contour_95)
}


#' Visualize Cell Neighbourhoods on Dimensionality Reduction
#'
#' @title visualize_neighbourhood
#' @description Visualizes the neighbours of cells belonging to a specific group
#'   on a dimensionality reduction plot (e.g., UMAP, t-SNE, PCA). Can display
#'   neighbours either as simple highlighted points or as density contours showing
#'   the spatial distribution patterns. This helps understand how cell types
#'   interact spatially and identify neighbourhood structures.
#' @param obj A Seurat, SingleCellExperiment or SCNeighbours object containing single-cell data.
#' If in Seurat or SingleCellExperiment form will firist be converted to SCneighbours format
#' @param meta_data_column Name of the metadata column in the object to pull values from
#'   for identifying the cells of interest.
#' @param meta_data_highlight The specific value within meta_data_column to
#'   analyze neighbours for (e.g., a specific cluster or cell type).
#' @param reduction Name of the dimensionality reduction to plot
#'   (e.g., "umap", "tsne", "pca").
#' @param density Logical indicating whether to plot density contours (TRUE) or
#'   simple point highlights (FALSE). Default is FALSE.
#' @param graph either a nearest neigbour graph in igraph, dgCMatrix or Seurat format, or the name of a graph stored in the Seurat object.
#'   (e.g., "RNA_nn", "RNA_snn", or "SCT_nn").
#' @param percent Percentile level for density contours when density=TRUE
#'   (default: 95). Only used when density=TRUE.
#' @return A ggplot2 object showing the dimensionality reduction with neighbours
#'   highlighted either as red points (density=FALSE) or as density contours
#'   (density=TRUE).
#' @export
#' @importFrom Seurat Embeddings
#' @importFrom dplyr mutate
#' @importFrom ggplot2 ggplot aes geom_point geom_path geom_density_2d theme_classic ggtitle
#' @importFrom rlang sym
#' @importFrom magrittr %>%
#' @importFrom Matrix colSums
visualize_neighbourhood = function(obj, meta_data_column, meta_data_highlight, reduction = NULL, density = F, graph = "RNA_nn", percent = 95) {
  # all neighbour cells
	scn <- check_single_cell_object(obj, graph, reduction)
	
  #n = colnames(seu@graphs[[graph]])[Matrix::colSums(seu@graphs[[graph]][seu[[meta_data_column]] == meta_data_highlight,]) > 0] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform
	
	g <- scn[['graph']]
	meta <- scn[['metadata']]
	emb <- scn[['embeddings']]
	axis <- scn[['key']]
	
	if(!meta_data_column %in% names(meta))
	{
		stop(paste(meta_data_column, " not found in metadata."))
	}
	
	n = colnames(g)[Matrix::colSums(g[meta[[meta_data_column]] == meta_data_highlight,]) > 0] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform
	
	
  # reduction emeddings of these neighbour cells
  # Embeddings(seu, reduction = reduction) %>%
  #   as_tibble(rownames = "bc") %>%
  #   mutate(neighbour = bc %in% n) -> d
	
	emb %>% mutate(neighbour = bc %in% n) -> d

  # matching axis names to the reduction names, please adjust it according to the real case
  #axis = emb@key
  #axis <- gsub("umap","UMAP",gsub("dm", "DC", reduction))
  axis_1 <- rlang::sym(paste0(axis, "1"))
  axis_2 <- rlang::sym(paste0(axis, "2"))

  if (!density){
    # plot by just highlighting these neighbour cells in the reduction map
    ggplot(d, aes(!!axis_1, !!axis_2)) +
      geom_point(colour = "grey") +
      geom_point(data = d[d$neighbour,], colour = "red", size = 0.6) +
      theme_classic()+ ggtitle(paste(meta_data_column, meta_data_highlight))

  } else{
    # d %>% mutate(m = seu[[meta_data_column]]) %>%
    #   group_by(m) %>% summarise(f = sum(neighbour)/n()) %>%
    #   filter(f>0.1) %>% as.data.frame() -> x
    # print(x)
    # print(as.list(levels(factor(x[,1]))))
    # plot the density plot of these neighbour cells in the reduction map
	  # contour_95 <- CalculateContour(seu, reduction = reduction,
	  #                                meta_data_column = meta_data_column, meta_data_highlight = meta_data_highlight,
	  #                                percent = percent)
  	
  	contour_95 <- CalculateContour(scn,
  	                                meta_data_column = meta_data_column, meta_data_highlight = meta_data_highlight,
  	                                percent = percent)

	  # allconts = lapply(as.list(c("CX3CR1 LN" , "TPEX LN" ,   "TPEX Tumor")), function(t)
	  #   CalculateContour(seu, meta_data_column = meta_data_column,
	  #                    meta_data_highlight = t, reduction = reduction,
	  #                    percent = 90) %>% mutate(g = t))
	  #
	  # allc = purrr::reduce(allconts, rbind)
	  # allc$g = factor(allc$g, levels = levels(seu_cd8_cdr3$CellType_major))

    ggplot(d, aes(!!axis_1, !!axis_2)) +
      geom_point(colour = "grey") +
      #geom_point(data = d[seu[[meta_data_column]] == meta_data_highlight,], size = 2, colour = "black")+
      #geom_point(data = d[seu[[meta_data_column]] == meta_data_highlight,], size = 1, colour = "grey")+
      geom_path(aes(x, y), data=contour_95, linetype = "dashed", linewidth = 1.5) +
      geom_density_2d(data = d[d$neighbour,]) +
      theme_classic() + ggtitle(paste0(meta_data_column, meta_data_highlight))
  }
}
