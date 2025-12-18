test_that("Seurat works", {
	expect_s3_class(SCneighbours::check_single_cell_object(seu, 'RNA_nn', 'umap.unintegrated'), 'scn_object')
	expect_s3_class(SCneighbours::check_single_cell_object(seu, sce.graph, 'umap.unintegrated'), 'scn_object')
})

test_that("SCE works", {
	expect_s3_class(SCneighbours::check_single_cell_object(sce, sce.graph, 'UMAP.UNINTEGRATED'), 'scn_object')
})


