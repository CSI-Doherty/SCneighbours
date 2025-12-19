test_that("Visulize neighbourhood for Seurat works", {
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 2, 'umap')))
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(seu, meta_data_column = 'seurat_clusters', meta_data_highlight = 2, 'umap', density = T)))
})

test_that("Visulize neighbourhood for SCE works", {
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(sce, meta_data_column = 'seurat_clusters', meta_data_highlight = 2, 'UMAP', graph = sce.graph)))
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(sce, meta_data_column = 'seurat_clusters', meta_data_highlight = 2, 'UMAP', graph = sce.graph, density = T)))
})

test_that("Visulize neighbourhood percentage for Seurat works", {
	expect_true(ggplot2::is_ggplot(visualise_neighbour_percentage(seu, meta_data_column = 'seurat_clusters',  graph = sce.graph)))
})

test_that("Visulize neighbourhood percentage for SCE works", {
	expect_true(ggplot2::is_ggplot(visualise_neighbour_percentage(sce, meta_data_column = 'seurat_clusters', graph = sce.graph)))
})

