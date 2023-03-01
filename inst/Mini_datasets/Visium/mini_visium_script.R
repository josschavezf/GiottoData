## MINI VISIUM script and dataset preparation ##


#remotes::install_github("drieslab/Giotto@suite_dev")
library(Giotto)

#remotes::install_github("drieslab/GiottoData")
library(GiottoData) # devtools::load_all()

library(data.table)

# 0. preparation ####
# ----------------- #

## create instructions
instrs = createGiottoInstructions(save_dir = tempdir(),
                                  save_plot = FALSE,
                                  show_plot = TRUE,
                                  return_plot = FALSE)

## provide path to vizgen folder
data_path = system.file('/Mini_datasets/Visium/Raw/', package = 'GiottoData')

## 0.1 path to images ####
# ---------------------- #
image_path = paste0(data_path, '/', 'images/deg_image.png')

## 0.2 path to expression matrix ####
# --------------------------- #
expr_path = paste0(data_path, '/', 'visium_DG_expr.txt.gz')

## 0.3 path to spot locations ####
# -------------------------------------- #
locations_path = paste0(data_path, '/', 'visium_DG_locs.txt')


# 1. create subcellular dataset with transcript and polygon information ####
# ------------------------------------------------------------------------ #

testmat = readExprMatrix(expr_path)

testmat = read_expression_data(expr_path)

seq_along(testmat)
testmat[[1]]


set_expression_values

testloc = read_spatial_location_data()

devtools::load_all('/Users/rubendries/Packages/R_Packages/Giotto')

mini_visium <- createGiottoObject(expression = expr_path,
                                  spatial_locs = locations_path,
                                  instructions = instrs)

showGiottoSpatLocs(mini_visium)
showGiottoExpression(mini_visium)

## 1.1. add image ####
# ------------------ #

spatlocsDT = getSpatialLocations(mini_visium)
mini_extent = terra::ext(c(range(spatlocsDT$sdimx), range(spatlocsDT$sdimy)))
imagelist = createGiottoLargeImageList(raster_objects = image_path,
                                       names = 'image',
                                       extent = mini_extent)
mini_visium = addGiottoImage(gobject = mini_visium,
                             largeImages = imagelist)
showGiottoImageNames(mini_visium)

## 1.2. visualize ####
# ------------------ #
spatPlot2D(gobject = mini_visium,
           spat_unit = 'cell',
           show_image = TRUE,
           largeImage_name = 'image',
           point_shape = 'no_border',
           point_size = 2.5,
           point_alpha = 0.4)


# 2 process ####
# ------------ #
mini_visium <- normalizeGiotto(gobject = mini_visium, scalefactor = 6000, verbose = T)

## filter
mini_visium <- filterGiotto(gobject = mini_visium,
                            expression_threshold = 1,
                            feat_det_in_min_cells = 5,
                            min_det_feats_per_cell = 20,
                            expression_values = c('raw'),
                            verbose = T)

## add gene & cell statistics
mini_visium <- addStatistics(gobject = mini_visium)

## visualize
spatPlot2D(gobject = mini_visium,
           show_image = TRUE,
           largeImage_name = 'image',
           point_alpha = 0.7)

spatPlot2D(gobject = mini_visium,
           show_image = TRUE,
           largeImage_name = 'image',
           point_alpha = 0.7,
           cell_color = 'nr_feats',
           color_as_factor = F)


# 3 dimension reduction ####
# ------------------------ #
mini_visium <- calculateHVF(gobject = mini_visium)

## run PCA on expression values (default)
mini_visium <- runPCA(gobject = mini_visium)
screePlot(mini_visium, ncp = 30)

plotPCA(gobject = mini_visium)

## run UMAP and tSNE on PCA space (default)
mini_visium <- runUMAP(mini_visium, dimensions_to_use = 1:10)
plotUMAP(gobject = mini_visium)

mini_visium <- runtSNE(mini_visium, dimensions_to_use = 1:10)
plotTSNE(gobject = mini_visium)

# 4 cluster ####
# ------------ #

## sNN network (default)
mini_visium <- createNearestNetwork(gobject = mini_visium, dimensions_to_use = 1:5, k = 10)
## Leiden clustering
mini_visium <- doLeidenCluster(gobject = mini_visium, resolution = 0.1, n_iterations = 1000)

plotUMAP(gobject = mini_visium, cell_color = 'leiden_clus', show_NN_network = T, point_size = 2.5)

spatDimPlot(gobject = mini_visium,
            show_image = TRUE,
            largeImage_name = 'image',
            cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)



## 9. save Giotto object ####
# ------------------------- #
format(object.size(mini_visium), units = 'Mb')

# you need to use your local GiottoData repo
giottodata_repo = '/Users/rubendries/Packages/R_Packages/GiottoData/inst/Mini_datasets/'

saveGiotto(mini_visium,
           foldername = 'VisiumObject',
           #dir = paste0(system.file(package = 'GiottoData'),'/', 'Mini_datasets/Vizgen/'),
           dir = paste0(giottodata_repo, '/', 'Visium/'),
           overwrite = TRUE)

pDataDT(mini_visium)
pDataDT(mini_visium)



## some quick tests ##
visium_test = loadGiotto(path_to_folder = system.file('/Mini_datasets/Visium/VisiumObject/',
                                                      package = 'GiottoData'))



spatDimPlot(gobject = visium_test,
            show_image = TRUE,
            largeImage_name = 'image',
            cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)






