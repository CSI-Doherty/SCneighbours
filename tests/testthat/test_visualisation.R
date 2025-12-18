test_that("Visulize neighbourhood for Seurat works", {
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(seu, meta_data_column = 'unintegrated_clusters', meta_data_highlight = 2, 'umap.unintegrated')))
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(seu, meta_data_column = 'unintegrated_clusters', meta_data_highlight = 2, 'umap.unintegrated', density = T)))
})




test_that("Visulize neighbourhood for SCE works", {
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(sce, meta_data_column = 'unintegrated_clusters', meta_data_highlight = 2, 'UMAP.UNINTEGRATED', graph = sce.graph)))
	expect_true(ggplot2::is_ggplot(visualize_neighbourhood(sce, meta_data_column = 'unintegrated_clusters', meta_data_highlight = 2, 'UMAP.UNINTEGRATED', graph = sce.graph, density = T)))
})
