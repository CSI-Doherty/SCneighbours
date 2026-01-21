test_that("Seurat neighbours distances are valid", {
	results <- SCneighbours::calculate_neighbour_distance_for_all_cells(seu, 'umap', graph = 'RNA_nn', colname = 'neighbour_distance')
	expect_type(results, 'Seurat')
	expect_equal(round(sum(results$neighbour_distance)), round(1151.0))
})

test_that("SCE neighbour distances are valid", {
	results <- SCneighbours::calculate_neighbour_distance_for_all_cells(sce, 'UMAP', graph = sce.graph, colname = 'neighbour_distance')
	expect_type(results, 'SingleCellExperiment')
	expect_equal(round(sum(results$neighbour_distance)), round(2888.0))
})


test_that("SCE and Seurat neighbour distances are equal", {
	results <- SCneighbours::calculate_neighbour_distance_for_all_cells(sce, 'UMAP', graph = seu.graph, colname = 'neighbour_distance')
	expect_type(results, 'SingleCellExperiment')
	expect_equal(round(sum(results$neighbour_distance)), round(1151.0))
})

