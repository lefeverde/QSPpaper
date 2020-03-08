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
