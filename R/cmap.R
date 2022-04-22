


#' Converts GCTX file to list signature list
#'
#' @param gctx_file gctx file
#' @param cids signatures to parse (see \code{\link{parse.gctx}})
#' @param decreasing_order whether to sort in decreasing order
#'
#' @importFrom cmapR parse_gctx
#'
#' @return signature list formatted for \code{\link[QSPpaper]{broad_cmap_score}}
#' @export
#'
#' @examples
gctx_to_sigList <- function(gctx_file, cids, decreasing_order=TRUE){
  sig_scores <- parse_gctx(gctx_file, cid = cids)
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






#' Calculates CMAP score using the Broad's latest methodology for a single instance
#'
#' This function calculates the CMAP score according to the method outlined  \code{\href{https://clue.io/connectopedia/cmap_algorithms}{here}} for a single instance of upregulated-genes, down-regulated genes and drug-signature
#'
#' @param up_genes vector of up-regulated genes
#' @param down_genes vector of down-regulated genes
#' @param drug_sig_vec drug signature vector
#'
#' @importFrom fgsea calcGseaStat
#'
#' @return
#' @export
#'
#' @examples
broad_cmap_score <- function(up_genes, down_genes, drug_sig_vec){

  # Checks if the vectors are sorted in Descending order by
  # sign fliping and checks if vals are from smallest to greatest
  assort_check <- !is.unsorted(drug_sig_vec*-1)
  stopifnot("drug_sig_vec is not in descending order" = assort_check)

  up_gene_idx <-
    match(up_genes, names(drug_sig_vec)) %>%
    na.omit
  down_gene_idx <-
    match(down_genes, names(drug_sig_vec)) %>%
    na.omit
  es_up <- calcGseaStat(drug_sig_vec, up_gene_idx)
  es_down <- calcGseaStat(drug_sig_vec, down_gene_idx)
  if(sign(es_up) == sign(es_down)){
    cmap_score <- 0
  } else {
    cmap_score <- (es_up - es_down)/2
  }
  out_df <- tibble(es_up, es_down, cmap_score)
  return(out_df)
}

#' Wrapper to run the \code{\link{broad_cmap_score}} for a single gene set
#'
#' This wrapper gets the CMAP scores for a single up and down-regulated gene set across  all drug signatures
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
single_broad_wrapper <- function(up_genes, down_genes, drug_sig_list, n_cores){



  es_list <- mclapply(seq_along(drug_sig_list), function(i){
    nm <- names(drug_sig_list)[i]
    x <- drug_sig_list[[i]]
    es_df <- broad_cmap_score(up_genes, down_genes, x)
    es_df <- cbind(sig_id=nm, es_df)
  },
  mc.cores = n_cores) %>%
    bind_rows
  up_genes <- paste(up_genes, collapse = ';')
  down_genes <- paste(down_genes, collapse = ';')
  es_list <- cbind(up_genes, down_genes, es_list)
  return(es_list)
}

.sort_check <- function(drug_sig_vec){
  out <- is.unsorted(rev(drug_sig_vec))
  return(out)
}

#' Wrapper to run the \code{\link{broad_cmap_score}} for many gene sets
#'
#' TODO Add more on the required input data format for gene signatures
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
multi_broad_wrapper <- function(query_genes, drug_sig_list, n_cores){
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
        single_broad_wrapper(up_genes, down_genes, drug_sig_list, n_cores) %>%
        cbind(nm, .)
    }) %>% do.call(rbind, .)
  return(cmap_gsea_res)
}





.cmap_build02_wrapper <- function(dz_genes_up, dz_genes_down, cmap_list, n_perm=1000, seed=NULL){
  #TODO fix signature inputs so an empty
  # signature can be passed

  cmap_score <-
    lapply(seq_along(cmap_list), function(i){
      cmap_exp_signature <-
        cmap_list[[i]]

      dz_cmap_scores <-
        cmap_score_new(dz_genes_up,
                       dz_genes_down,
                       cmap_exp_signature)
    }) %>% do.call(c, .)

  if(!is.null(seed)){
    set.seed(seed)
  }
  # Random selection of drug signatures
  random_sig_ids <-
    sample(seq_along(cmap_list), n_perm, replace=TRUE)
  # genes in cmap data
  gene_ids <- cmap_list[[1]]$ids
  dz_up_num <-
    dz_genes_up[dz_genes_up %in% gene_ids] %>%
    na.omit %>%
    length

  dz_down_num <-
    dz_genes_down[dz_genes_down %in% gene_ids] %>%
    na.omit %>%
    length

  dz_total_num <- dz_up_num + dz_down_num


  random_cmap_res <-
    lapply(random_sig_ids, function(idx){

      cmap_exp_signature <-
        cmap_list[[idx]]


      random_input_signature_genes <-
        sample(gene_ids, dz_total_num, replace=FALSE)


      if(dz_up_num == 0){
        rand_dz_gene_up <- c(NA)
      } else {
        rand_dz_gene_up <-
          random_input_signature_genes[1:dz_up_num]
      }

      if(dz_down_num == 0){
        rand_dz_gene_down <- c(NA)
      } else {
        rand_dz_gene_down <-
          random_input_signature_genes[!random_input_signature_genes %in% rand_dz_gene_up]
        # Changed from index because it got confusing
      }

      dz_cmap_scores <-
        cmap_score_new(rand_dz_gene_up,
                       rand_dz_gene_down,
                       cmap_exp_signature)
    }) %>% do.call(c, .)

  p_values <-
    lapply(cmap_score, function(score){
      sum(abs(random_cmap_res) >= abs(score))/length(random_cmap_res)
    }) %>% do.call(c, .)

  fdr_p_value <- p.adjust(p_values, method = 'fdr')


  cmap_res <-
    data.frame(pert_id=names(cmap_list),
               cmap_score,
               p_values,
               fdr_p_value)

  return(cmap_res)
}


#' Internal function to create a data.frame with the ordered permuations
#'
#' @inheritParams broad_cmap_score
#'
#' @param n_perm number of permutations to perform
#' @param seed random number seed
#' @param num_chunks number of chunk levels for chunk_idx column
#'
#'
#' @importFrom cmapR read_gctx_ids
#' @return
#' @export
#'
#' @examples
make_permute_sort_df <- function(cids, n_perm=1000, num_chunks=4, seed=NULL){
  if(!is.null(seed)){
    set.seed(seed)
  }
  sig_df <-
    cids %>%
    enframe %>%
    setNames(., c("idx", "sig_id"))
  perm_ids <-
    sample(sig_df$sig_id, size = n_perm, replace = TRUE)
  id_cnts <- table(perm_ids) %>%
    as.matrix %>%
    data.frame(., check.rows = FALSE, check.names = FALSE) %>%
    rownames_to_column %>%
    setNames(., c("sig_id", "count"))
  out_df <-
    inner_join(sig_df, id_cnts, by="sig_id") %>%
    arrange(idx)
  chunk_idx <-
    seq(0, nrow(out_df) - 1) %% num_chunks %>%
    sort
  out_df <- cbind(chunk_idx, out_df)

  return(out_df)
}

#' Runs the permutation testing to get a vector of CMAP scores
#'
#' @inheritParams gctx_to_sigList
#' @inheritParams broad_cmap_score
#' @inheritParams single_broad_wrapper
#'
#' @importFrom cmapR read_gctx_meta
#' @return
#' @export
#'
#' @examples
get_random_cmap_scores <- function(up_genes, down_genes, gctx_file, cids=NULL, n_perm=1000, num_chunks=4, seed=NULL, n_cores=3){
  if(!is.null(seed)){
    set.seed(seed)
  }
  gene_ids <-
    read_gctx_ids(gctx_file, "row")
  if(is.null(cids)){
    cids <- read_gctx_ids(gctx_file, "col")
  }

  comb_genes <- c(up_genes, down_genes)
  if(!all(comb_genes %in% gene_ids)){
    stop("Not all gene signature genes are in the pert signatures.")
  }
  sorted_permute_df <-
    make_permute_sort_df(cids, n_perm = n_perm, num_chunks = num_chunks, seed = seed) %>%
    split(., .$chunk_idx)

  cmap_scores <- map(sorted_permute_df, function(x){
    # loads chunk
    sig_list <-
      gctx_to_sigList(gctx_file, x$sig_id, decreasing_order = TRUE)
    # Use sig name for key
    sig_ids <- names(sig_list)
    random_scores <- mclapply(sig_ids, function(nm){
      pert_vec <- sig_list[[nm]]
      n <- x %>%
        filter(sig_id == nm) %>%
        pull(count)
      # Runs the CMAP score function with a random gene set
      # according to the number determined by make_permute_sort_df
      tmp_scores <- map(1:n, function(i){
        random_genes <-
          sample(gene_ids, length(comb_genes), replace=FALSE)
        r_up_genes <-
          random_genes[1:length(up_genes)]
        r_down_genes <-
          random_genes[!random_genes %in% r_up_genes]
        cs <- broad_cmap_score(r_up_genes, r_down_genes, pert_vec)
        # Adds in the random gene signature not sure if this is
        # a good idea or not.
        r_up_genes <- paste(r_up_genes, collapse = ';')
        r_down_genes <- paste(r_down_genes, collapse = ';')
        cs <- cbind(r_up_genes, r_down_genes, cs)
      }) %>% bind_rows
    }, mc.cores = n_cores) %>% bind_rows
  }) %>% bind_rows
  # Binds everything into a single df
  return(cmap_scores)

}


#' Calculates p-values for CMap results
#'
#' This function requires the output CMap results as created by \code{\link{single_broad_wrapper}} and a
#' matrix of permuation values created by \code{\link{get_random_cmap_scores}}. The output is a data.frame
#' of the CMap results with the p-values added.
#'
#'
#'
#' @param cmap_res output results from \code{\link{single_broad_wrapper}}
#' @param perm_res output results from \code{\link{get_random_cmap_scores}}
#'
#' @return
#' @export
#'
#' @examples
add_p_vals_to_cmap_single <- function(cmap_res, perm_res){
  res_w_p_vals <- seq_along(cmap_res) %>%
    map(function(idx){
      x <- cmap_res[[idx]]
      y <- perm_res[[idx]]
      # Ok so the function outer actually made things slower than map
      # so I'll reduce the computations per iter. That might help
      # but I'm not expecting much. I'm also parallelizing it which
      # might help.
      x_scores <- abs(x$cmap_score) # res cmap scores
      y_scores <- abs(y$cmap_score) # random cmap scores
      y_length <- length(y$cmap_score)
      p_value <- x_scores %>%
        map_dbl(function(score){
          sum(y_scores >= score)/y_length
        })
      fdr_p_value <- p.adjust(p_value, method = 'fdr')
      out_df <- cbind(x, p_value, fdr_p_value)
    }) %>% bind_rows
  return(res_w_p_vals)
}



#' Same as \code{\link{add_p_vals_to_cmap_single}} but parallel
#'
#' @inheritParams add_p_vals_to_cmap_single
#'
#' @return
#' @export
#'
#' @examples
add_p_vals_to_cmap_para <- function(cmap_res, perm_res){
  res_w_p_vals <- seq_along(cmap_res) %>%
    map(function(idx){
      x <- cmap_res[[idx]]
      y <- perm_res[[idx]]
      # Ok so the function outer actually made things slower than map
      # so I'll reduce the computations per iter. That might help
      # but I'm not expecting much. I'm also parallelizing it which
      # might help.
      x_scores <- abs(x$cmap_score) # res cmap scores
      y_scores <- abs(y$cmap_score) # random cmap scores
      y_length <- length(y$cmap_score)
      p_value <- x_scores %>%
        future_map_dbl(function(score){
          sum(y_scores >= score)/y_length
        })
      fdr_p_value <- p.adjust(p_value, method = 'fdr')
      out_df <- cbind(x, p_value, fdr_p_value)
    }) %>% bind_rows
  return(res_w_p_vals)
}


#' Normalizes the connectivity score by cell type
#'
#' Scores are normalized by cell type.
#'
#'
#' @param cs_mat data.frame of CS scores
#'
#' @return
#' @export
#'
#' @examples
calculate_ncs <- function(cs_mat){
  # There's not enough ncs sigs for all the different cell types so
  # I'm not just not going to bother with it for the time being.


  mu_pos_df <- cs_mat %>%
    filter(cmap_score > 0) %>%
    # filter(is_ncs_sig == 1) %>%
    # filter(is_touchstone == 1) %>%
    group_by(sig_idx, cell_iname) %>%
    summarise(mu_pos=mean(cmap_score)) %>%
    ungroup


  mu_neg_df <- cs_mat %>%
    filter(cmap_score < 0) %>%
    # filter(is_ncs_sig == 1) %>%
    # filter(is_touchstone == 1) %>%
    group_by(sig_idx, cell_iname) %>%
    summarise(mu_neg=mean(cmap_score)) %>%
    ungroup

  out_df <- left_join(cs_mat, mu_pos_df) %>% left_join(., mu_neg_df)
  out_df <- out_df %>%
    mutate(norm_cmap_score=if_else(cmap_score > 0, cmap_score/mu_pos, -cmap_score/mu_neg)) %>%
    select(-mu_neg, -mu_pos)
  # This is done to set the divide by 0 values to 0
  out_df$norm_cmap_score[is.na(out_df$norm_cmap_score)] <- 0
  return(out_df)
}

#' Calculates pert score per cell type
#'
#' @inheritParams calculate_ncs
#' @param quant_thresh percentile threshold. Default is .67 which was used in Subramanian 2017
#' @param quant_type This is how the quantiles are calculated. See \link[stats]{quantile}
#'
#' @return
#' @export
#'
#' @examples
calculate_pert_cell_score <- function(cs_mat, quant_thresh=.67, quant_type=5){
  hi_thresh <- quant_thresh
  lo_thresh <- 1 - quant_thresh
  # Doign a group_by with 3 vars seems to slow the code way down
  # splitting is faster
  out_df <- split(cs_mat, cs_mat$sig_idx) %>%
    map_dfr(function(x){
      sig_idx <- unique(x$sig_idx)
      pert_cell_mat <- x %>%
        group_by(pert_id, cell_iname) %>%
        summarise(pert_cell_hi=quantile(norm_cmap_score, probs=hi_thresh, type=quant_type),
                  pert_cell_lo=quantile(norm_cmap_score, probs=lo_thresh, type=quant_type),
                  n_sigs=n()) %>%
        ungroup %>%
        mutate(pert_cell_sum_score=if_else(abs(pert_cell_hi) > abs(pert_cell_lo), pert_cell_hi, pert_cell_lo))
      pert_cell_mat <- cbind(sig_idx, pert_cell_mat)
    })

  return(out_df)

}


#' Calculates the perturbagen summary score
#'
#' @param pert_cell_mat output from \link{calculate_pert_cell_score}
#' @inheritParams calculate_pert_cell_score
#'
#' @return
#' @export
#'
#' @examples
calculate_pert_score <- function(pert_cell_mat, quant_thresh=.67, quant_type=5){
  hi_thresh <- quant_thresh
  lo_thresh <- 1 - quant_thresh
  # Doign a group_by with 3 vars seems to slow the code way down
  # splitting is faster
  out_df <- split(pert_cell_mat, pert_cell_mat$sig_idx) %>%
    map_dfr(function(x){
      sig_idx <- unique(x$sig_idx)
      pert_mat <- x %>%
        group_by(pert_id) %>%
        summarise(pert_hi=quantile(pert_cell_sum_score, probs=hi_thresh, type=quant_type),
                  pert_lo=quantile(pert_cell_sum_score, probs=lo_thresh, type=quant_type),
                  n_sigs=n()) %>%
        ungroup %>%
        mutate(pert_sum_score=if_else(abs(pert_hi) > abs(pert_lo), pert_hi, pert_lo))
      pert_mat <- cbind(sig_idx, pert_mat)
    })
  out_df <- out_df %>%
    group_by(sig_idx) %>%
    arrange(pert_sum_score, .by_group = TRUE) %>%
    mutate(drug_idx=row_number()) %>%
    ungroup

  return(out_df)

}


#' Returns the top 20 compounds using the best score approach
#'
#' @param best_cs_mat output from \link{best_score_mat}
#'
#' @return
#' @export
#'
#' @examples
best_score_top20 <- function(best_cs_mat){
  top20_mat <-  best_cs_mat %>%
    distinct %>%
    group_by(sig_idx) %>%
    arrange(cmap_score, .by_group=TRUE) %>%
    slice(1:20) %>%
    ungroup
  top20_mat <- top20_mat %>%
    group_by(pert_iname) %>%
    summarise(freq=n(), unique_inst=length(unique(sig_id))) %>%
    arrange(desc(freq))

  return(top20_mat)
}

#' Chooses the best score per compound
#'
#' @inheritParams calculate_ncs
#'
#' @return
#' @export
#'
#' @examples
best_score_mat <- function(cs_mat){
  best_cs_mat <- cs_mat %>%
    filter(cmap_score < 0) %>%
    filter(fdr_p_value < .05) %>%
    group_by(sig_idx, pert_iname) %>%
    arrange(cmap_score, .by_group=TRUE) %>%
    slice(1) %>%
    ungroup %>%
    group_by(sig_idx) %>%
    arrange(cmap_score, .by_group=TRUE) %>%
    mutate(drug_idx=row_number()) %>%
    ungroup
}

#' Returns the top 20 compounds after performing maximum quantile statistic
#'
#' @param pert_mat
#'
#' @return
#' @export
#'
#' @examples
lincs2017_score_top20 <- function(pert_mat){
  clean_mat <- pert_mat %>%
    group_by(sig_idx, pert_iname) %>%
    arrange(pert_sum_score, .by_group = TRUE) %>%
    slice(1) %>%
    ungroup

  top20_mat <-  clean_mat %>%
    distinct %>%
    group_by(sig_idx) %>%
    arrange(pert_sum_score, .by_group=TRUE) %>%
    slice(1:20, with_ties=FALSE) %>%
    ungroup %>%
    group_by(pert_iname) %>%
    summarise(freq=n()) %>%
    arrange(desc(freq))
  return(top20_mat)

}
