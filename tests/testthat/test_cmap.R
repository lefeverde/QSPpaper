library(QSPpaper)
library(testthat)
context('Testing functions from cmap.R')

test_gctx <-
  system.file('extdata',
              'test_sig_500x978.gctx_n500x978.gctx',
              package = 'QSPpaper',
              mustWork = TRUE)
test_data <-
  system.file('extdata',
              'test_gene_set.RData',
              package = 'QSPpaper',
              mustWork = TRUE)
test_res <-
  system.file('extdata',
              'test_cmap_res.rds',
              package = 'QSPpaper',
              mustWork = TRUE)

load(test_data)
cids <- cmapR::read_gctx_ids(test_gctx, "col")
drug_vec_list <- gctx_to_sigList(test_gctx, cids)

test_res <- readRDS(test_res)

test_sig <-
  seq(-10, 10) %>%
  sort(., decreasing = FALSE) %>%
  setNames(., paste0("G", seq(0, 20)))

test_that("broad_cmap_score works when up & down gene sets give same score ", {
  up_genes <- c("G0","G1","G2")
  down_genes <- c("G20", "G19", "G18")
  e_res <- t(c(1, -1, 1)) %>%
    as_tibble(., .name_repair = "unique") %>%
    setNames(., c("es_up", "es_down","cmap_score"))
  o_res <- broad_cmap_score(up_genes, down_genes, test_sig)
  expect_equal(as.data.frame(e_res), as.data.frame(o_res), tolerance=.001)
})


test_that("broad_cmap_score works when up has higher rank ", {
  up_genes <- c("G0","G1","G2")
  down_genes <- c("G12", "G13", "G14")
  e_res <- t(c(1, -0.667, 0.833)) %>%
    as_tibble(., .name_repair = "unique") %>%
    setNames(., c("es_up", "es_down","cmap_score"))
  o_res <- broad_cmap_score(up_genes, down_genes, test_sig)
  expect_equal(as.data.frame(e_res), as.data.frame(o_res), tolerance=.001)
})

test_that("broad_cmap_score works when down has higher rank ", {
  up_genes <- c("G6", "G7", "G8")
  down_genes <- c("G20", "G19", "G18")
  e_res <- t(c(0.667, -1, 0.833)) %>%
    as_tibble(., .name_repair = "unique") %>%
    setNames(., c("es_up", "es_down","cmap_score"))
  o_res <- broad_cmap_score(up_genes, down_genes, test_sig)
  expect_equal(as.data.frame(e_res), as.data.frame(o_res), tolerance=.001)
})

test_that("broad_cmap_score returns 0 when up and down yeild the same sign ", {
  up_genes <- c("G0", "G19", "G16", "G15"  )
  down_genes <- c("G20","G18", "G17", "G1", "G2")
  e_res <- t(c(-0.490, -0.533, 0)) %>%
    as_tibble(., .name_repair = "unique") %>%
    setNames(., c("es_up", "es_down","cmap_score"))
  o_res <- broad_cmap_score(up_genes, down_genes, test_sig)
  expect_equal(as.data.frame(e_res), as.data.frame(o_res), tolerance=.001)
})

test_that("single_broad_wrapper returns the correct CMap results",{
  o_res <-
    single_broad_wrapper(up_genes, down_genes, drug_vec_list, 1)
  o_res <- o_res %>% arrange(cmap_score)
  e_res <- test_res %>% arrange(cmap_score)
  expect_equal(as.data.frame(e_res), as.data.frame(o_res), tolerance=.001)

})
