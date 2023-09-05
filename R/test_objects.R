
#' @title Load a mini giotto subobject
#' @name loadSubObjectMini
#' @description Utility function to load a mini giotto suite S4 subobject.
#' Can be useful for testing functions quickly or providing contained examples.
#' @param x subobject type to load
#' @param idx which of multiple example subobjects to load when more than one
#' is available (see \code{\link{listSubObjectMini}})
#' @export
loadSubObjectMini = function(x, idx = 1L) {

  avail_obj_dt = list_subobject_mini()
  data_path = avail_obj_dt[type == x & index == idx, path]

  load_data = readRDS(file = data_path)

  # wrapped objects
  if(x %in% c('giottoPoints', 'giottoPolygon')) {
    load_data = GiottoClass::vect(load_data)
  }

  if(x == 'giottoLargeImage') {
    original_path = load_data@file_path
    new_path = gsub(pattern = '.*[/]GiottoData/', replacement = '', x = original_path)
    new_path = paste0(gDataDir(), new_path)
    load_data = GiottoClass::reconnect_giottoLargeImage(giottoLargeImage = load_data,
                                                        image_path = new_path)
  }

  if(x == 'giottoImage') {
    original_path = load_data@file_path
    new_path = gsub(pattern = '.*[/]GiottoData/', replacement = '', x = original_path)
    new_path = paste0(gDataDir(), new_path)
    load_data = GiottoClass::reconnect_giottoImage_MG(giottoImage = load_data,
                                                      image_path = new_path)
  }

  return(load_data)

}


#' @title List available mini giotto subobjects
#' @name listSubObjectMini
#' @description Lists the available mini giotto subobjects to load. Items are first
#' classified by the type of subobject, then by the index number which is used
#' to tell apart files when there is more than one example subobject of that type
#' available. The final column is the file name.
#' @param x subobject type (NULL lists all subobject types)
#' @export
listSubObjectMini = function(x = NULL) {

  avail_obj_dt = list_subobject_mini()
  avail_obj_dt_show = copy(avail_obj_dt)[, path := NULL]
  data.table::setcolorder(avail_obj_dt_show, neworder = c('type', 'index', 'file'))

  avail_obj_dt_show = data.table::copy(avail_obj_dt_show)

  if(is.null(x)) return(avail_obj_dt_show)
  else return((avail_obj_dt_show[type == x,]))

}


#' @title List available subobjects data
#' @name list_subobject_mini
#' @description Utility function to list the available saved subobjects. This is
#' the internal function which returns a version of the table with the full pathfiles
#' and is not prettied or subsetted by type. See \code{\link{listSubObjectMini}}
#' for the external version which provides those functions.
#' @keywords internal
list_subobject_mini = function() {

  miniobj_path = paste0(gDataDir(), '/Mini_objects/subobjects')

  avail_obj = list.files(path = miniobj_path, full.names = TRUE, recursive = TRUE)
  avail_obj_dt = data.table::data.table(file = basename(avail_obj), path = avail_obj)
  avail_obj_dt[, 'type' := gsub(pattern = paste0(miniobj_path, '/'), replacement = '', x = path)]
  avail_obj_dt[, 'type' := gsub(pattern = '[/].*', replacement = '', x = type)]

  avail_obj_dt[, 'index' := seq_along(path), by = type]

  return(avail_obj_dt)

}




