#' @title loadGiottoMini
#' @name loadGiottoMini
#' @param dataset mini dataset giotto object to load
#' @param python_path pythan path to use
#' @description This function will automatically load one of the existing mini
#' giotto objects. These giotto objects can be used to test Giotto functions
#' and run examples. If no python path is provided it will try to find and use
#' the Giotto python environment. Images associated with the giotto mini objects
#' will be reconnected if possible.
#'
#' Instructions, such as for saving plots, can be changed
#' using the \code{\link{changeGiottoInstructions}}
#' @export
loadGiottoMini = function(dataset = c('visium', 'seqfish', 'starmap', 'vizgen'),
                          python_path = NULL) {


  dataset = match.arg(dataset, choices = c('visium', 'seqfish', 'starmap', 'vizgen'))


  if(dataset == 'visium') {

    # 1. load giotto object
    mini_gobject = readRDS(system.file("/Mini_datasets/Visium/gobject_mini_visium.RDS", package = 'GiottoData'))

    # 2. add image back
    image_path = system.file("/Mini_datasets/Visium/images/deg_image.png", package = 'Giotto')
    spatlocsDT = Giotto::get_spatial_locations(mini_gobject)
    mini_extent = terra::ext(c(range(spatlocsDT$sdimx), range(spatlocsDT$sdimy)))
    imagelist = Giotto::createGiottoLargeImageList(raster_objects = image_path,
                                                   names = 'image',
                                                   extent = mini_extent)
    mini_gobject = Giotto::addGiottoImage(gobject = mini_gobject,
                                          largeImages = imagelist)

  }


  if(dataset == 'vizgen') {

    mini_gobject = Giotto::loadGiotto(path_to_folder = system.file('/Mini_datasets/Vizgen/VizgenObject/', package = 'GiottoData'),
                              python_path = python_path)

  }


  if(dataset == 'seqfish') {
    Giotto:::wrap_msg('To be implemented \n')
  }

  if(dataset == 'starmap') {
    Giotto:::wrap_msg('To be implemented \n')
  }


  # 1. change default instructions
  identified_python_path = Giotto:::set_giotto_python_path(python_path = python_path)
  mini_gobject = Giotto::changeGiottoInstructions(gobject = mini_gobject,
                                          params = c('python_path', 'show_plot', 'return_plot', 'save_plot', 'save_dir'),
                                          new_values = c(identified_python_path, TRUE, FALSE, FALSE, NA))

  return(mini_gobject)

}


#' @title getSpatialDataset
#' @name getSpatialDataset
#' @param dataset dataset to download
#' @param directory directory to save the data to
#' @param verbose verbosity
#' @param dryrun dryrun: does not download data but shows download commands
#' @param \dots additional parameters to \code{\link[utils]{download.file}}
#' @description This function will automatically download the spatial locations and
#' expression matrix for the chosen dataset. These files are already in the right format
#' to create a Giotto object. If wget is installed on your machine, you can add
#' 'method = wget' to the parameters to download files faster.
#' @export
getSpatialDataset = function(dataset = c('ST_OB1',
                                         'ST_OB2',
                                         'codex_spleen',
                                         'cycif_PDAC',
                                         'starmap_3D_cortex',
                                         'osmfish_SS_cortex',
                                         'merfish_preoptic',
                                         'mini_seqFISH',
                                         'seqfish_SS_cortex',
                                         'seqfish_OB',
                                         'slideseq_cerebellum',
                                         'ST_SCC',
                                         'scRNA_prostate',
                                         'scRNA_mouse_brain'),
                             directory = getwd(),
                             verbose = TRUE,
                             dryrun = FALSE,
                             ...) {

  sel_dataset = match.arg(arg = dataset, choices = c('ST_OB1',
                                                     'ST_OB2',
                                                     'codex_spleen',
                                                     'cycif_PDAC',
                                                     'starmap_3D_cortex',
                                                     'osmfish_SS_cortex',
                                                     'merfish_preoptic',
                                                     'mini_seqFISH',
                                                     'seqfish_SS_cortex',
                                                     'seqfish_OB',
                                                     'slideseq_cerebellum',
                                                     'ST_SCC',
                                                     'scRNA_prostate',
                                                     'scRNA_mouse_brain'))

  # check operating system first
  os_specific_system = Giotto:::get_os()

  #if(os_specific_system == 'windows') {
  #  stop('This function is currently not supported on windows systems,
  #       please visit https://github.com/RubD/spatial-datasets and manually download your files')
  #}


  # check directory
  if(!file.exists(directory)) {
    warning('The output directory does not exist and will be created \n')
    dir.create(directory, recursive = T)
  }

  datasets_file = system.file("extdata", "datasets.txt", package = 'GiottoData')
  datasets_file = data.table::fread(datasets_file, sep = "\t")



  ## check if wget is installed
  #message = system("if ! command -v wget &> /dev/null
  #                  then
  #                  echo 'wget could not be found, please install wget first'
  #                  exit
  #                  fi", intern = TRUE)

  #if(identical(message, character(0))) {
  #  print('wget was found, start downloading datasets: ')
  #} else {
  #  stop(message)
  #}

  ## alternative
  #wget_works = try(system('command -v wget', intern = T))

  #if(class(wget_works) == 'try-error' | is.na(wget_works[1])) {
  #  stop('wget was not found, please install wget first \n')
  #} else {
  #  print('wget was found, start downloading datasets: \n')
  #}


  selection = datasets_file[['dataset']] == sel_dataset
  selected_dataset_info = datasets_file[selection,]
  if(verbose) {
    Giotto:::wrap_msg('Selected dataset links for: ', sel_dataset, ' \n \n')
    print(selected_dataset_info)
  }



  # get url to expression matrix and download
  if(verbose) {
    Giotto:::wrap_msg("\n \n Download expression matrix: \n")
  }
  expr_matrix_url = selected_dataset_info[['expr_matrix']]

  if(expr_matrix_url == "") {
    Giotto:::wrap_msg('\n No expression found, skip this step \n')
  } else {

    expr_matrix_url = unlist(strsplit(expr_matrix_url, split = '\\|'))

    for(url in expr_matrix_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)

      if(dryrun) {
        Giotto:::wrap_msg("utils::download.file(url = ", url, ", destfile = ", mydestfile, ", ...)")
      } else {
        utils::download.file(url = url, destfile = mydestfile, ...)
      }
    }

  }






  # get url to spatial locations and download
  if(verbose) {
    Giotto:::wrap_msg("\n \n Download spatial locations: \n")
  }

  spatial_locs_url = selected_dataset_info[['spatial_locs']]

  if(spatial_locs_url == "") {
    Giotto:::wrap_msg('\n No spatial locations found, skip this step \n')
  } else {

    spatial_locs_url = unlist(strsplit(spatial_locs_url, split = '\\|'))

    for(url in expr_matrix_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)

      if(dryrun) {
        Giotto:::wrap_msg("utils::download.file(url = ", url, ", destfile = ", mydestfile, ", ...)")
      } else {
        utils::download.file(url = url, destfile = mydestfile, ...)
      }
    }
  }





  # get url(s) to additional metadata files and download
  if(verbose) {
    Giotto:::wrap_msg("\n \n Download metadata: \n")
  }

  #metadata_url = selected_dataset_info[['metadata']][[1]]
  metadata_url = selected_dataset_info[['metadata']]

  if(metadata_url == "") {
    Giotto:::wrap_msg('\n No metadata found, skip this step \n')
  } else {

    metadata_url = unlist(strsplit(metadata_url, split = '\\|'))

    for(url in metadata_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)

      if(dryrun) {
        Giotto:::wrap_msg("utils::download.file(url = ", url, ", destfile = ", mydestfile, ", ...)")
      } else {
        utils::download.file(url = url, destfile = mydestfile, ...)
      }
    }

  }

}
