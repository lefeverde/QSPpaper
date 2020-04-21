#' Creates a data.frame to make a GSEA plot
#'
#' This function creates a data.frame which can be used to
#' create the line plot of the enrichment score. It was taken
#' and slightly modified from the \code{\link[fgsea]{plotEnrichment}}.
#'
#' @param query_idx indexes of query which correpond to position in \code{\link{ordered_vec}}
#' @param ordered_vec ordered vector of scores
#' @param gsea_weight power to raise the ordered vector (see \code{\href{https://www.genepattern.org/modules/docs/GSEAPreranked/1}{here}} for details.)
#'
#' @importFrom fgsea calcGseaStat
#'
#' @return
#' @export
#'
#' @examples
create_gsea_plot_data <- function(query_idx, ordered_vec, gsea_weight=0){
  query_idx <- sort(query_idx)
  # Spaghet because of wierd way R does exponents
  ordered_vec <-
    sign(ordered_vec)*(abs(ordered_vec)**gsea_weight)
  es <-
    calcGseaStat(ordered_vec, query_idx, returnAllExtremes = TRUE)

  n <- length(ordered_vec)
  xs <- as.vector(rbind(query_idx - 1, query_idx))
  ys <- as.vector(rbind(es$bottoms, es$tops))
  toPlot <-
    data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  return(toPlot)
}


#' Creates a running sum GSEA plot
#'
#' This is my own customized version of the plot created
#' by \code{\link[fgsea]{plotEnrichment}}.
#'
#' @param plot_data data.frame created by \code{\link{create_gsea_plot_data}}
#' @param main_title optional title to be added to plot
#'
#' @importFrom cowplot theme_cowplot insert_xaxis_grob
#'
#' @return
#' @export
#'
#' @examples
make_gsea_plot <- function(plot_data, main_title=""){

  y_max <-
    max(abs(plot_data$y))
  x_max <-
    plot_data$x[which(abs(plot_data$y) == y_max)]

  lbl <- paste0("Max ES: ", formatC(y_max, 2))

  line_plt <-
    ggplot(data=plot_data, aes(x, y)) +
    geom_line() +
    theme_cowplot() +
    scale_y_continuous(limits = c(-1, 1)) +
    scale_x_continuous(position = "bottom") +
    geom_hline(yintercept = 0, linetype="dashed", colour="blue") +
    geom_segment(x=x_max, xend=x_max, y=0, yend=y_max, linetype="dotted", colour="red") +
    geom_text(aes(label=lbl, x=x_max, y=y_max), nudge_y=.1) +
    theme(axis.text = element_text(face="bold"),
          axis.title = element_text(face="bold"),
          plot.title = element_text(hjust = .5))

  line_plt <- line_plt +
    labs(x="Rank ordered index", y="Enrichment score (ES)",title = main_title)

  tick_idxs <- plot_data[2:(nrow(plot_data) -1),] %>%
    mutate(plot_tick=rep_len(c(TRUE, FALSE), nrow(.)))

  tick_plt <- tick_idxs %>%
    filter(plot_tick) %>%
    ggplot(data=.) +
    geom_vline(aes(xintercept = x)) +
    theme_void()

  tmp_plt <-
    insert_xaxis_grob(line_plt,
                      tick_plt,
                      position = "top",
                      height = grid::unit(0.1, "null"))
  return(tmp_plt)
}

#' Creates legend paramaters for \code{\link[ComplexHeatmap]{HeatmapAnnotation}}
#'
#' The reason I made a function here is that I could never remember how the
#' parameter list is supposed to be formatted. So it's somewhat superfluous except
#' that this is hard to google.
#'
#' @param tit title
#' @param leg_size legend size (in mm)
#' @param leg_width legend width (in mm)
#' @param leg_height legend height  (in mm)
#' @param font_size font size
#' @param direction direction of legend (either horizontal or vertical )
#'
#' @importFrom grid unit gpar
#'
#' @return
#' @export
#'
#' @examples
create_annot_leg_params <- function(tit=NULL, leg_size=7.5, leg_width=NULL, leg_height=NULL, font_size=14, direction='vertical'){
  if(is.null(leg_width)){
    leg_width <- leg_size
  }
  if(is.null(leg_height)){
    leg_height <- leg_size
  }
  leg_params <- list(
    grid_height = unit(leg_height, "mm"),
    grid_width = unit(leg_width, "mm"),
    labels_gp = gpar(fontsize = font_size),
    title_gp = gpar(fontsize = font_size, fontface='bold'),
    legend_direction=direction)
  if(!is.null(tit)){
    leg_params[['title']] <- tit
  }
  return(leg_params)

}

#' Creates legend paramaters for \code{\link[ComplexHeatmap]{Heatmap}}
#'
#' Same sort of function as \code{\link[QSPpaper]{create_annot_leg_params}}
#'
#' @inheritParams create_annot_leg_params
#'
#' @return
#' @export
#'
#' @examples
create_heatmap_leg_params <- function(tit=NULL, leg_size=50, font_size=14, direction='horizontal'){
  leg_params <- list(

    legend_width = unit(leg_size, "mm"),
    legend_height = unit(leg_size, "mm"),
    labels_gp = gpar(fontsize = font_size),
    title_gp = gpar(fontsize = font_size, fontface='bold'),
    # title_position = 'leftcenter-rot',
    legend_direction=direction)
  if(!is.null(tit)){
    leg_params[['title']] <- tit
  }
  return(leg_params)

}
