test_that("SCE outside neighbours are valid", {
	results <- SCneighbours::calculate_outside_neighbours_cell(sce, 'seurat_clusters', graph = sce.graph, colname = 'outside')
	expect_s4_class(results, 'SingleCellExperiment')
	expect_equal(sum(results$outside), 84055.169)
})

test_that("Seurat outside neighbours are valid", {
	results <- SCneighbours::calculate_outside_neighbours_cell(seu, 'seurat_clusters', colname = 'outside')
	expect_s4_class(results, 'Seurat')
	expect_equal(sum(results$outside), 32180)
})

test_that("SCE and Seurat outside neighbours are equal", {
	expect_equal(sum(SCneighbours::calculate_outside_neighbours_cell(sce, 'seurat_clusters', graph = seu.graph, colname = 'outside')$outside), 32180)
})

test_that("Seurat neighbour_percentage_all_ids are valid", {
  expect_s3_class(SCneighbours::calculate_neighbour_percentage_all_ids(seu, 'seurat_clusters'), "data.frame")
})


test_that("SCE neighbour_percentage_all_ids are valid", {
  expect_s3_class(SCneighbours::calculate_neighbour_percentage_all_ids(sce, 'seurat_clusters', graph = sce.graph), "data.frame")
})

test_that("Seurat calculate_neighbour_percentage are valid", {
  expect_s3_class(SCneighbours::calculate_neighbour_percentage(seu,  meta_data_column = 'seurat_clusters', meta_data_highlight = 2), "data.frame")
})

test_that("SCE calculate_neighbour_percentage are valid", {
  expect_s3_class(SCneighbours::calculate_neighbour_percentage(seu,  meta_data_column = 'seurat_clusters', meta_data_highlight = 2, graph = sce.graph), "data.frame")
})