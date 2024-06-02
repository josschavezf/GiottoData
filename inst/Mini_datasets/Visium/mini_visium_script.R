## MINI VISIUM script and dataset preparation ##

library(Giotto)
library(GiottoData)

# NOTE:
# deconvolution and HMRF do not work too well on a small subset of the full
# dataset.
# Deconvolution expects all cell types to be present when starting
# from a particular reference single cell dataset. In a subset, this is going
# to be less true, producing strange and noisy results.
# HMRF just looks neater when run on a full tissue.

# load in a pre-made object based on the visium mouse brain tutorial
g <- readRDS("???")
# replace images downstream...



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
image_path <- vector("list")
image_path["alignment"] = file.path(data_path, 'images/deg_image.png')
image_path["image"] =  file.path(data_path, 'images/deg_image2.png')

## 0.2 path to expression matrix ####
# --------------------------- #
expr_path = paste0(data_path, '/', 'visium_DG_expr.txt.gz')

## 0.3 path to spot locations ####
# -------------------------------------- #
locations_path = paste0(data_path, '/', 'visium_DG_locs.txt')

## 0.4 path to metadata ####
# -------------------------------------- #
meta_path = paste0(data_path, '/', 'visium_DG_meta.txt')

## 0.5 path to scalefactors ####
# -------------------------------------- #
scalef_path <- file.path(data_path, "scalefactors_json.json")


# 1. create visium dataset ####
# ------------------------------------------------------------------------ #
mini_visium <- createGiottoObject(expression = expr_path,
                                  spatial_locs = locations_path,
                                  cell_metadata = meta_path,
                                  instructions = instrs)

mini_visium <- addVisiumPolygons(mini_visium, scalefactor_path = scalef_path)
g <- addVisiumPolygons(g, scalefactor_path = scalef_path)

showGiottoSpatLocs(mini_visium)
showGiottoExpression(mini_visium)

## 1.1. add image ####
# ------------------ #

spatlocsDT = getSpatialLocations(mini_visium)
mini_extent = terra::ext(c(range(spatlocsDT$sdimx), range(spatlocsDT$sdimy)))
image_align = createGiottoLargeImage(
  raster_object = image_path$alignment,
  name = "alignment",
  extent = terra::ext(2364.5, 6522.5, -5425.25, -2620.75)
 )
image_he <- createGiottoLargeImage(
    raster_object = image_path$image,
    name = "image",
    extent = terra::ext(2000.5, 6790.5, -5730.25, -2380.75)
)
# image_he <- rescale(image_he, 13.2420382165605, 13.2287735849057)
# image_he <- spatShift(image_he, dx = 4322, dy = -3937)
imagelist <- list(image_align, image_he)
names(imagelist) <- c("alignment", "he")
mini_visium = addGiottoImage(gobject = mini_visium,
                             images = imagelist)
g <- addGiottoImage(g, images = imagelist)
showGiottoImageNames(mini_visium)

# save mini
giottodata_repo <- GiottoData:::gdata_dataset_devdir()
saveGiotto(g,
           foldername = 'VisiumObject',
           dir = paste0(giottodata_repo, '/', 'Visium/'),
           overwrite = TRUE)


## 1.2. visualize ####
# ------------------ #
spatPlot2D(gobject = mini_visium,
           spat_unit = 'cell',
           show_image = TRUE,
           image_name = 'alignment',
           point_shape = 'no_border',
           point_size = 2.5,
           point_alpha = 0.4)

spatPlot2D(gobject = mini_visium,
           spat_unit = 'cell',
           show_image = TRUE,
           image_name = 'image',
           point_shape = 'no_border',
           point_size = 2.5,
           point_alpha = 0.4)

spatInSituPlotPoints(
  mini_visium,
  show_polygon = TRUE,
  polygon_color = "cyan",
  show_image = TRUE,
  image_name = "alignment"
)

mini_visium@images$alignment@raster_object

# 2 process ####
# ------------ #
mini_visium <- normalizeGiotto(gobject = mini_visium, scalefactor = 6000, verbose = T)

list_expression(mini_visium)
list_spatial_locations(mini_visium)

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
           image_name = 'image',
           point_alpha = 0.7)

spatPlot2D(gobject = mini_visium,
           show_image = TRUE,
           image_name = 'image',
           background_color = "black",
           cell_color_gradient = c("cyan", "blue", "black", "orange", "yellow"),
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

plotUMAP(gobject = mini_visium,
         cell_color = 'leiden_clus',
         show_NN_network = T,
         point_size = 2.5)

spatDimPlot(gobject = mini_visium,
            show_image = TRUE,
            image_name = 'image',
            cell_color = 'leiden_clus',
            dim_point_size = 2,
            spat_point_size = 2.5)


# 5. spatial network ####
# --------------------- #
mini_visium <- createSpatialNetwork(gobject = mini_visium)

mini_visium <- createSpatialNetwork(gobject = mini_visium,
                                    method = 'kNN', k = 10,
                                    maximum_distance_knn = 400,
                                    name = 'spatial_network')


# 6. spatial genes ####
# ------------------- #
showGiottoSpatNetworks(mini_visium)
ranktest <- binSpect(
    mini_visium,
    bin_method = 'rank',
    calc_hub = T,
    hub_min_int = 5,
    spatial_network_name = "Delaunay_network"
)


# 7. spatial co-expression ####
# --------------------------- #

# 7.1 cluster the top 500 spatial genes into 20 clusters
ext_spatial_genes = ranktest[1:300,]$feats

# here we use existing detectSpatialCorGenes function to calculate pairwise distances between genes (but set network_smoothing=0 to use default clustering)
spat_cor_netw_DT = detectSpatialCorFeats(mini_visium,
                                         method = 'network',
                                         spatial_network_name = 'spatial_network',
                                         subset_feats = ext_spatial_genes)

# 7.2 identify most similar spatially correlated genes for one gene
top10_genes = showSpatialCorFeats(spat_cor_netw_DT, feats = 'Dsp', show_top_feats = 10)

spatFeatPlot2D(mini_visium,
               expression_values = 'scaled',
               feats = top10_genes$variable[1:4], point_size = 3)


# 7.3 identify potenial spatial co-expression
spat_cor_netw_DT = clusterSpatialCorFeats(spat_cor_netw_DT, name = 'spat_netw_clus', k = 7)

# visualize clusters
heatmSpatialCorFeats(mini_visium,
                     spatCorObject = spat_cor_netw_DT,
                     use_clus_name = 'spat_netw_clus',
                     heatmap_legend_param = list(title = NULL),
                     save_param = list(base_height = 6, base_width = 8, units = 'cm'))


# 7.4 create metagenes / co-expression modules
cluster_genes = getBalancedSpatCoexpressionFeats(spat_cor_netw_DT, maximum = 30)
mini_visium = createMetafeats(mini_visium, feat_clusters = cluster_genes, name = 'cluster_metagene')

spatCellPlot(mini_visium,
             spat_enr_names = 'cluster_metagene',
             cell_annotation_values = as.character(c(1:7)),
             point_size = 1, cow_n_col = 3)




# 8. spatially informed clusters ####
# --------------------------------- #
my_spatial_genes = names(cluster_genes)

mini_visium <- runPCA(gobject = mini_visium,
                      feats_to_use = my_spatial_genes,
                      name = 'custom_pca')
mini_visium <- runUMAP(mini_visium,
                       dim_reduction_name = 'custom_pca',
                       dimensions_to_use = 1:20,
                       name = 'custom_umap')
mini_visium <- createNearestNetwork(gobject = mini_visium,
                                    dim_reduction_name = 'custom_pca',
                                    dimensions_to_use = 1:20, k = 5,
                                    name = 'custom_NN')
mini_visium <- doLeidenCluster(gobject = mini_visium,
                               network_name = 'custom_NN',
                               resolution = 0.15, n_iterations = 1000,
                               name = 'custom_leiden')

spatPlot2D(mini_visium,
           cell_color = 'custom_leiden')




# 9. DWLS #
# ------- #

# DWLS does not run well on a small dataset without the full representation
# of cells from a mouse brain. Results will be merged in from pre-made
# Refer to the Visium Mouse Brain tutorial on the Giotto website for
# how the DWLS spatial enrichment object was generated


# 10. save Giotto object ####
# ------------------------- #
format(object.size(mini_visium), units = 'Mb')

# you need to use your local GiottoData repo
#giottodata_repo = '/Users/rubendries/Packages/R_Packages/GiottoData/inst/Mini_datasets/'
#giottodata_repo = '/Users/rubendries/r_packages/GiottoData//inst/Mini_datasets/'

saveGiotto(mini_visium,
           foldername = 'VisiumObject',
           #dir = paste0(system.file(package = 'GiottoData'),'/', 'Mini_datasets/Vizgen/'),
           dir = paste0(giottodata_repo, '/', 'Visium/'),
           overwrite = TRUE)

pDataDT(mini_visium)


## some quick tests ##
visium_test = loadGiotto(path_to_folder = system.file('/Mini_datasets/Visium/VisiumObject/',
                                                      package = 'GiottoData'))


spatPlot2D(visium_test,
           show_image = T,
           image_name = 'image',
           cell_color = 'custom_leiden')

spatDimPlot(gobject = visium_test,
            show_image = TRUE,
            image_name = 'image',
            cell_color = 'leiden_clus',
            dim_point_size = 2,
            spat_point_size = 2.5)




