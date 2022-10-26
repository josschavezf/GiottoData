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
    mini_gobject = readRDS(system.file("/Mini_datasets/Visium/gobject_mini_visium.RDS", package = 'Giotto'))

    # 2. add image back
    image_path = system.file("/Mini_datasets/Visium/images/deg_image.png", package = 'Giotto')
    spatlocsDT = Giotto::get_spatial_locations(mini_gobject)
    mini_extent = terra::ext(c(range(spatlocsDT$sdimx), range(spatlocsDT$sdimy)))
    imagelist = createGiottoLargeImageList(raster_objects = image_path,
                                                   names = 'image',
                                                   extent = mini_extent)
    mini_gobject = addGiottoImage(gobject = mini_gobject,
                                          largeImages = imagelist)

  }


  if(dataset == 'vizgen') {

    # 1. load giotto object
    mini_gobject = readRDS(system.file("/Mini_datasets/Vizgen/gobject_mini_vizgen.RDS", package = 'Giotto'))

    # 2. add spatvectors back in place (giottoPoints and giottoPolygons)
    rna_spatVector = terra::vect(system.file("/Mini_datasets/Vizgen/processed_data/rna_spatVector.shp", package = 'Giotto'))
    mini_gobject@feat_info$rna@spatVector = rna_spatVector

    z0_spatVector = terra::vect(system.file("/Mini_datasets/Vizgen/processed_data/z0_spatVector.shp", package = 'Giotto'))
    mini_gobject@spatial_info$z0@spatVector = z0_spatVector
    z0_spatVectorCentroids = terra::vect(system.file("/Mini_datasets/Vizgen/processed_data/z0_spatVectorCentroids.shp", package = 'Giotto'))
    mini_gobject@spatial_info$z0@spatVectorCentroids = z0_spatVectorCentroids

    z1_spatVector = terra::vect(system.file("/Mini_datasets/Vizgen/processed_data/z1_spatVector.shp", package = 'Giotto'))
    mini_gobject@spatial_info$z1@spatVector = z1_spatVector
    z1_spatVectorCentroids = terra::vect(system.file("/Mini_datasets/Vizgen/processed_data/z1_spatVectorCentroids.shp", package = 'Giotto'))
    mini_gobject@spatial_info$z1@spatVectorCentroids = z1_spatVectorCentroids



    # 3. add spatRaster back in place (largeGiottoImages)

    # x and y information from original script
    ultra_mini_extent = terra::ext(c(6400.029, 6900.037, -5150.007, -4699.967 ))

    # location
    DAPI_z0_image_path = system.file("/Mini_datasets/Vizgen/images/mini_dataset_dapi_z0.jpg", package = 'Giotto')
    DAPI_z1_image_path = system.file("/Mini_datasets/Vizgen/images/mini_dataset_dapi_z1.jpg", package = 'Giotto')
    polyT_z0_image_path = system.file("/Mini_datasets/Vizgen/images/mini_dataset_polyT_z0.jpg", package = 'Giotto')
    polyT_z1_image_path = system.file("/Mini_datasets/Vizgen/images/mini_dataset_polyT_z1.jpg", package = 'Giotto')

    # create image list
    image_paths = c(DAPI_z0_image_path, DAPI_z1_image_path,
                    polyT_z0_image_path, polyT_z1_image_path)
    image_names = c('dapi_z0', 'dapi_z1',
                    'polyT_z0', 'polyT_z1')

    imagelist = createGiottoLargeImageList(raster_objects = image_paths,
                                           names = image_names,
                                           negative_y = TRUE,
                                           extent = ultra_mini_extent)

    mini_gobject = addGiottoImage(gobject = mini_gobject,
                                  largeImages = imagelist)

  }


  if(dataset == 'seqfish') {
    wrap_msg('To be implemented \n')
  }

  if(dataset == 'starmap') {
    wrap_msg('To be implemented \n')
  }


  # 1. change default instructions
  identified_python_path = set_giotto_python_path(python_path = python_path)
  mini_gobject = changeGiottoInstructions(gobject = mini_gobject,
                                          params = c('python_path', 'show_plot', 'return_plot', 'save_plot', 'save_dir'),
                                          new_values = c(identified_python_path, TRUE, FALSE, FALSE, NA))

  return(mini_gobject)

}


#' @title getSpatialDataset
#' @name getSpatialDataset
#' @param dataset dataset to download
#' @param directory directory to save the data to
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
                                         'Human_PCa_scRNAseq',
                                         'Mouse_brain_scRNAseq'),
                             directory = getwd(),
                             ...) {

  sel_dataset = match.arg(dataset, choices = c('ST_OB1',
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
                                               'Human_PCa_scRNAseq',
                                               'Mouse_brain_scRNAseq'))

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
  datasets_file = data.table::fread(datasets_file)



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



  # get url to spatial locations and download
  spatial_locs_url = datasets_file[dataset == sel_dataset][['spatial_locs']]
  myfilename = basename(spatial_locs_url)
  mydestfile = paste0(directory,'/', myfilename)

  print(spatial_locs_url)
  print(mydestfile)

  utils::download.file(url = spatial_locs_url, destfile = mydestfile, ...)

  #system(paste0("wget -P ", "'",directory,"'"," ", spatial_locs_url))


  # get url to expression matrix and download
  expr_matrix_url = datasets_file[dataset == sel_dataset][['expr_matrix']]
  myfilename = basename(expr_matrix_url)
  mydestfile = paste0(directory,'/', myfilename)
  utils::download.file(url = expr_matrix_url, destfile = mydestfile, ...)

  #system(paste0("wget -P ", "'",directory,"'"," ", expr_matrix_url))

  # get url(s) to additional metadata files and download
  metadata_url = datasets_file[dataset == sel_dataset][['metadata']][[1]]
  metadata_url = unlist(strsplit(metadata_url, split = '\\|'))

  if(identical(metadata_url, character(0))) {
    NULL
  } else {
    for(url in metadata_url) {
      myfilename = basename(url)
      mydestfile = paste0(directory,'/', myfilename)
      utils::download.file(url = url, destfile = mydestfile, ...)
      #system(paste0("wget -P ", "'",directory,"'"," ", url))
    }
  }

}
