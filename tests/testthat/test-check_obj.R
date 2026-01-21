test_that("Seurat works", {
	expect_s3_class(SCneighbours::check_single_cell_object(seu, 'RNA_nn', 'umap'), 'scn_object')
	expect_s3_class(SCneighbours::check_single_cell_object(seu, sce.graph, 'umap'), 'scn_object')
	expect_s3_class(SCneighbours::check_single_cell_object(seu), 'scn_object')
})

test_that("SCE works", {
	expect_s3_class(SCneighbours::check_single_cell_object(sce, sce.graph, 'UMAP'), 'scn_object')
	expect_s3_class(SCneighbours::check_single_cell_object(sce, sce.graph), 'scn_object')
})

test_that("Seurat graphs are valid", {
  expect_type(SCneighbours::check_single_cell_object(seu, sce.graph, 'umap')[['graph']], "dgCMatrix")
	expect_type(SCneighbours::check_single_cell_object(seu, seu.graph)[['graph']], "dgCMatrix")
	expect_type(SCneighbours::check_single_cell_object(seu)[['graph']], "dgCMatrix")
})


test_that("SCE graphs are valid", {
	expect_type(SCneighbours::check_single_cell_object(sce, seu.graph, 'UMAP')[['graph']], "dgCMatrix")
	expect_type(SCneighbours::check_single_cell_object(sce, sce.graph, 'UMAP')[['graph']], "dgCMatrix")
})
