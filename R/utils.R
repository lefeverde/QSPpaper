#' Splits a vector into chunks
#'
#' @param x vector
#' @param num_chunks number of chunks
#'
#' @return list of chunked vectors
#' @export
#'
#' @examples
chunker <- function(x, num_chunks=3){
  chunks <-  split(x, sort((seq_along(x) - 1)%%num_chunks))
  return(chunks)
}

#' Uniquefys a res df by group by taking absolute max
#'
#' @param res_df a tidy results data.frame (see \code{\link[basicOmics]{get_limma_results}})
#' @param group_col  column with duplicates
#' @param val_col column of unique values
#'
#' @importFrom basicOmics uniquefy_by_abs_max
#'
#' @return a uniquefyd results data.frame
#' @export
#'
#' @examples
uniquefy_res_by_group <- function(res_df, group_col='external_gene_name', val_col='GeneRank'){
  tmp_res <- na.omit(res_df)
  tmp_res$GeneRank <- tmp_res$logFC* -log10(tmp_res$adj.P.Val)
  tmp_res <- tmp_res %>%
    group_by(., coefficient) %>%
    do(.,
       uniquefy_by_abs_max(., group_col, val_col))
  return(tmp_res)
}
