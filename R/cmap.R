


#' Converts GCTX file to list signature list
#'
#' @param gctx_file gctx file
#' @param cids signatures to parse (see \code{\link{parse.gctx}})
#' @param decreasing_order whether to sort in decreasing order
#'
#' @importFrom cmapR parse.gctx
#'
#' @return signature list formatted for \code{\link[QSPpaper]{broad_cmap_score}}
#' @export
#'
#' @examples
gctx_to_sigList <- function(gctx_file, cids, decreasing_order=TRUE){
  sig_scores <- parse.gctx(gctx_file, cid = cids)
  rns <- row.names(sig_scores@mat)
  sig_scores <- as_tibble(sig_scores@mat)


  sig_list <-
    as.list(sig_scores) %>%
    lapply(., function(x){
      x <- x %>%
        setNames(., rns)
      x <- sort(x, decreasing = decreasing_order)
    }) %>%
    setNames(., colnames(sig_scores))
  return(sig_list)
}




#' Calculates CMAP score according to \href{https://clue.io/connectopedia/cmap_algorithms}{CLUE}
#'
#' @param up_genes vector of up-regulated genes
#' @param down_genes vector of down-regulated genes
#' @param drug_sig_list list of drug signatures
#' @param n_cores number of cores
#'
#' @importFrom fgsea calcGseaStat
#' @importFrom parallel mclapply
#'
#'
#'
#' @return
#' @export
#'
#' @examples
broad_cmap_score <- function(up_genes, down_genes, drug_sig_list, n_cores){
  es_list <- mclapply(drug_sig_list, function(x){
    up_gene_idx <-
      match(up_genes, names(x)) %>%
      na.omit
    down_gene_idx <-
      match(down_genes, names(x)) %>%
      na.omit
    es_up <- calcGseaStat(x, up_gene_idx)
    es_down <- calcGseaStat(x, down_gene_idx)
    if(sign(es_up) == sign(es_down)){
      es <- 0
    } else {
      es <- (es_up - es_down)/2
    }
    return(es)
  },
  mc.cores = n_cores)
  out_df <-
    tibble(sig_id=names(es_list),
           cmap_score=unlist(es_list))
  return(out_df)
}


#' Wrapper to run the \code{\link{broad_cmap_score}} for many gene sets
#'
#' @param query_genes data.frame of gene signatures
#' @inheritParams broad_cmap_score
#' @param n_cores number of cores
#' @importFrom stringr str_split
#'
#'
#' @return
#' @export
#'
#' @examples
cmap_fgsea_wrapper <- function(query_genes, drug_sig_list, n_cores){
  cmap_gsea_res <-
    lapply(1:nrow(query_genes), function(i){
      nm <- query_genes[i,1]
      up_genes <-
        query_genes[i,2] %>%
        str_split(., ';') %>%
        unlist
      down_genes <-
        query_genes[i,3] %>%
        str_split(., ';') %>%
        unlist
      tmp_res <-
        broad_cmap_score(up_genes, down_genes, drug_sig_list, n_cores) %>%
        cbind(nm, .)
    }) %>% do.call(rbind, .)
  return(cmap_gsea_res)
}
