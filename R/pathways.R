#' Wrapper to do GSEA on tidy results
#'
#' @param res_df a tidy results data.frame of genes (see \code{\link[basicOmics]{get_limma_results}})
#' @param nPerm number of permutations to run
#' @param pvalueCutoff p-value cutoff of output to retain
#'
#' @importFrom  clusterProfiler gseKEGG
#'
#' @return a tidy results data.frame of KEGG pathways
#' @export
#'
#' @examples
kegg_gsea_by_coef <- function(res_df, nPerm=10000, pvalueCutoff=.1){

  sres <- uniquefy_res_by_group(res_df, 'entrezgene') %>%
    split(., .$coefficient)

  ol <- lapply(sres, function(x){
    print(unique(x$coefficient))
    gene_rank <- x$t %>%
      setNames(., x$entrezgene) %>%
      sort(., decreasing = TRUE)
    gsea_res <-
      gseKEGG(gene_rank,
              keyType = 'ncbi-geneid',
              nPerm = nPerm,
              pvalueCutoff = pvalueCutoff)
    return(gsea_res)
  })
  return(ol)
}
