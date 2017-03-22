#' c.AAbin
#' @export
c.AAbin <- function(..., recursive = FALSE)
{
  if (!all(unlist(lapply(list(...), is.list))))
    stop("the 'c' method for \"AAbin\" accepts only lists")
  structure(NextMethod("c"), class = "AAbin")
}
