#' Wrapper to create a fit object (see \code{\link[limma]{eBayes}}) using the contrast method
#'
#' @param v a voom object
#' @param group_cont vector of contrasts
#' @param mod model matrix
#'
#' @return fit object
#' @export
#'
#' @examples
make_cont_fit <- function(v, group_cont, mod){
  m <- data.frame(mod)
  cont_mod <- makeContrasts(contrasts = group_cont, levels=m)
  fit <- lmFit(v, m) %>%
    contrasts.fit(., cont_mod) %>%
    eBayes(.)
  return(fit)
}
