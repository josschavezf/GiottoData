#' @title Get GiottoData library path
#' @name gDataDir
#' @description Utility function to find library install location of the GiottoData
#' package. Provided as shorthand to make loading files easier.
#' @export
gDataDir <- function() {
    system.file(package = "GiottoData")
}
