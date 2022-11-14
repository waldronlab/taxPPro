
#' Choose column for attribute values
#'
#' \code{chooseColVal} chooses the column where attribute values are
#' actually stored.
#'
#' @param df A data frame imported with bugphyzz.
#'
#' @return A character string.
#' @export
#'
chooseColVal <- function(df) {

    av <- df$Attribute_value

    if(is.logical(av)) {
        attr_col <- 'Attribute'
    } else if (is.numeric(av) || is.character(av)) {
        attr_col <- 'Attribute_value'
    }

    attr_col
}
