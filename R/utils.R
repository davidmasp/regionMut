## Utils and helpers




#' Unlist with same name
#'
#' It's similar to unlist but here it takes the names and puts them in the
#' resulting vector as they are without appending numbers.
#'
#' @param x a list
#'
#' @return
#' # tI think this shouldn't be exported
#' @export
#'
#' @examples
regionmut_unlist <- function(x) {
  recode_vec = unlist(x,use.names = F)
  new_levels_rle = list(
    lengths = lengths(x),
    values = names(x)
  )
  names(recode_vec) = inverse.rle(new_levels_rle)
  return(recode_vec)
}

