#library(Seurat)
#library(dplyr)
#library(tidyverse)
#
## calculate the variance in coordinates of the neighbour cells for aech cell in a certain dimensionality reduction
#
##'@title neighbour_distance.scaled.
##'@description The function to calculate the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
##'@param i The cell index.
##'@param reduction The reduction map used to calculate the coordinate variance.
##'@param seu The seurat object.
##'@return A number, the average variance in coordinates of neighbour cells of a certain cell i based on the provided reduction map.
#neighbour_distance.scaled = function(i, reduction, seu, graph = "RNA_nn") {
#  # extract neighbour cells of cell i
#  n = colnames(seu@graphs[[graph]])[seu@graphs[[graph]][i,] == 1] # please change "RNA_nn" to "SCT_nn" if you are using SCTranform
#
#  # normalize dimension scales
#  embed.scale <- Embeddings(seu, reduction = reduction)
#  embed.scale[,1] <- BBmisc::normalize(embed.scale[,1], method = "range", range = c(-8, 8))
#  embed.scale[,2] <- BBmisc::normalize(embed.scale[,2], method = "range", range = c(-8, 8))
#
#  # dims of 20 cells
#  e = embed.scale[n, ]
#
#  # variance in dim
#  mean(var(e[,1]), var(e[,2]))
#}
#
#
#
##'@title calculate_neighbour_distance_for_all_cells.
##'@description The function to calculate the average variance in coordinates of neighbour cells of all cells in the seurat object based on the provided reduction map.
##'@param colname The colname to store the neighbourhood distance value in metadata.
##'@param reduction The reduction map used to calculate the coordinate variance.
##'@param seu The seurat object.
##'@return A seurat object with a new metadata column storing the varaince in coodinates of neighbour cells for each cell
#calculate_neighbour_distance_for_all_cells <- function(seu, reduction, colname, graph) {
#  for(i in 1:nrow(seu@meta.data)) {seu@meta.data[[colname]][i] = neighbour_distance.scaled(i, reduction, seu, graph)}
#  seu
#}
#
#calculate_neighbour_percentage <- function(seu, meta_data_column, meta_data_highlight, graph){
#	sn = colnames(seu@graphs[[graph]])[colSums(seu@graphs[[graph]][seu[[meta_data_column]] == meta_data_highlight,]) > 0]
#	ids = factor(seu@meta.data[sn,meta_data_column], levels = levels(factor(seu@meta.data[,meta_data_column])))
#	table(ids) %>% as.data.frame() %>%
#		mutate(f = Freq/sum(Freq)*100)
#}
#
#
#calculate_neighbour_percentage_all_ids <- function(seu, meta_data_column, graph){
#	ids = levels(factor(seu@meta.data[,meta_data_column]))
#	results = data.frame(ids = ids)
#	for(i in ids){
#		results[,i] = calculate_neighbour_percentage(seu, meta_data_column = meta_data_column, meta_data_highlight = i, graph)$f
#	}
#	return(results)
#}
#
#visualise_neighbour_percentage <- function(seu, meta_data_column, graph) {
#	x = calculate_neighbour_percentage_all_ids(seu, meta_data_column, graph)
#	d = dist(t(x[,-1]))
#	h = hclust(d)
#	
#	x %>% pivot_longer(-ids) %>% 
#		mutate(name = factor(name, levels = h$labels[h$order]),
#					 ids = factor(ids, levels = h$labels[h$order])) %>%
#		ggplot(aes(name,ids)) +
#		geom_tile(aes(fill = value)) +
#		scale_fill_gradientn(colours =  plotto_color_gradient_blue_red) +
#		theme_classic() +
#		theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#		labs(x = element_blank(), y=element_blank(), fill = "% shared")
#}
#
#calculate_outside_neighbours_cell <- function(seu, meta_data_column, graph, colname){
#	for(i in 1:nrow(seu@meta.data)){
#		ids = seu@meta.data[[meta_data_column]][seu@graphs[[graph]][i,]>0]
#		seu@meta.data[[colname]][i] = (1-sum(ids == seu@meta.data[[meta_data_column]][i])/length(ids))*100
#	}
#	return(seu)
#}
#


