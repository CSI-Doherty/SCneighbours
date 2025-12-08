library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyverse)

# utility function for contour plotting

CalculateContour = function(seu, meta_data_column, meta_data_highlight, reduction = "umap", percent = 95) {
  d = Embeddings(seu, reduction = reduction) %>% as_tibble(rownames = "bc")

  kd <- ks::kde(d[seu[[meta_data_column]] == meta_data_highlight,2:3], compute.cont=TRUE)
  contour_95 <- with(kd, contourLines(x=eval.points[[1]], y=eval.points[[2]],
                                      z=estimate, levels=cont[paste0(100-percent,"%")])[[1]])
  contour_95 <- data.frame(contour_95)
  #
  # d[seu[[meta_data_column]] == meta_data_highlight,2:3]
  return(contour_95)
}

# plot neighbour cells of an individual cell or a certain cell group (i.e., cluster)

#'@title visualize_neighbourhood.
#'@description The function to highlight the neighbour cells of certain cells.
#'@param meta_data_column Name of the column in seurat_object@meta.data slot to pull value from for highlighting.
#'@param meta_data_highlight Name of variable(s) within meta_data_name to highlight in the plot.
#'@param reduction The reduction map used to calculate the coordinate variance.
#'@param seu The seurat object.
#'@param density Whether to plot the density countour plot of neighbour cells. Default is FALSE.
#'@return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
visualize_neighbourhood = function(seu, meta_data_column, meta_data_highlight, reduction, density = F, graph = "RNA_nn", percent = 95) {
  # all neighbour cells
  n = colnames(seu@graphs[[graph]])[colSums(seu@graphs[[graph]][seu[[meta_data_column]] == meta_data_highlight,]) > 0] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform

  # reduction emeddings of these neighbour cells
  Embeddings(seu, reduction = reduction) %>%
    as_tibble(rownames = "bc") %>%
    mutate(neighbour = bc %in% n) -> d

  # matching axis names to the reduction names, please adjust it according to the real case
  axis = seu@reductions[[reduction]]@key
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
	  contour_95 <- CalculateContour(seu, reduction = reduction,
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
      geom_path(aes(x, y), data=contour_95, linetype = "dashed", size = 1.5) +
      geom_density_2d(data = d[d$neighbour,]) +
      theme_classic() + ggtitle(paste0(meta_data_column, meta_data_highlight))
  }
}

