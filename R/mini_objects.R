
#' @name generate_mini_subobjects
#' @title Generate mini subobjects
#' @description
#' Generate the GiottoData mini subobjects based on Mini_objects_script.R
#' in `/inst/Mini_objects/` subdirectory.
#' @returns NULL invisibly
#' @export
generate_mini_subobjects <- function() {

    script_path <- file.path(
        gdata_devdir(), "inst", "Mini_objects", "Mini_objects_script.R"
    )

    source(script_path, echo = TRUE)
    return(invisible())
}
