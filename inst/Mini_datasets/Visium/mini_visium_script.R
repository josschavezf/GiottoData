## MINI VISIUM script and dataset preparation ##

# devtools::load_all('/Users/rubendries/r_packages/Giotto/')
# devtools::load_all()

#remotes::install_github("drieslab/Giotto@suite_dev")
library(Giotto)

#remotes::install_github("drieslab/GiottoData")
library(GiottoData) # devtools::load_all()

library(terra)
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


# 9. save Giotto object ####
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



# 10. build from scratch ####
# -------------------------- #

devtools::load_all('/Users/rubendries/r_packages/Giotto/')

# 10.1 get expression data
list_expression(mini_visium)
raw_matrix = getExpression(mini_visium, values = 'raw', output = 'matrix')
normalized_matrix = getExpression(mini_visium, values = 'normalized', output = 'matrix')
scaled_matrix = getExpression(mini_visium, values = 'scaled', output = 'matrix')

# 10.2 get spatial location data
list_spatial_locations(mini_visium)
spatial_locations = getSpatialLocations(mini_visium, name = 'raw', output = 'data.table')

# 10.3 get cell and feature metadata
list_cell_metadata(mini_visium)
cell_metadata = getCellMetadata(mini_visium, output = 'data.table')
list_feat_metadata(mini_visium)
feat_metadata = getFeatureMetadata(mini_visium, output = 'data.table')

# 10.4 dimension reduction
list_dim_reductions(mini_visium)
pca_dim = getDimReduction(mini_visium, reduction_method = 'pca', name = 'pca', output = 'matrix')
pca_dim = getDimReduction(mini_visium, reduction_method = 'umap', name = 'umap', output = 'matrix')
pca_dim = getDimReduction(mini_visium, reduction_method = 'tsne', name = 'tsne', output = 'matrix')

# 10.5 get nearest networks
list_nearest_networks(mini_visium)
sNN_network = getNearestNetwork(mini_visium, nn_type = 'sNN', name = 'sNN.pca', output = 'data.table')

# 10.6 large images
list_images(mini_visium)
images = getGiottoImage(mini_visium, image_type = 'largeImage', name = 'image')
plot(images)


## 10.1 create spatialExperiment ####

speg = giottoToSpatialExperiment(giottoObj = mini_visium)

class(speg[[1]])

speg_rna = speg[[1]]

SummarizedExperiment::assayNames(speg_rna)

speg_rna@NAMES

speg_rna@assays@data$normalized_rna_cell

guard_against_notgiotto

?SpatialExperiment::SpatialExperiment()

# remake
mini_visium_remake <- createGiottoObject(expression = list('cell' =
                                                             list('rna' =
                                                                    list('raw' = raw_matrix,
                                                                         'normalized' = normalized_matrix,
                                                                         'scaled' = scaled_matrix))),
                                         instructions = instrs)


mini_visium_remake = setSpatialLocations(gobject = mini_visium_remake, spatlocs = spatial_locations)

mini_visium_remake = set_cell_metadata(gobject = mini_visium_remake,
                                       metadata = cell_metadata,
                                       spat_unit = 'cell', feat_type = 'rna')

mini_visium_remake = set_feature_metadata(gobject = mini_visium_remake,
                                       metadata = feat_metadata,
                                       spat_unit = 'cell', feat_type = 'rna')


read_dimension_reduction()

create_dim_obj()

mini_visium_remake = setDimReduction(gobject = mini_visium_remake, dimObject = )

spatDimPlot(gobject = mini_visium_remake,
            show_image = TRUE,
            largeImage_name = 'image',
            cell_color = 'leiden_clus',
            dim_point_size = 2, spat_point_size = 2.5)

str(mini_visium, 2)




