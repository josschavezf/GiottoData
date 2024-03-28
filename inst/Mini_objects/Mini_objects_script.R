## MINI OBJECTS script for use in subobject testing ##

# remotes::install_github("drieslab/GiottoData")
library(GiottoData) # devtools::load_all()

# remotes::install_github("drieslab/Giotto@suite_dev")

# 0. preparation ####
# ----------------- #

viz <- loadGiottoMini("vizgen")





# 1. extract subobjects ####
# ----------------- #

spat_locs <- getSpatialLocations(
    gobject = viz,
    spat_unit = "aggregate",
    name = "raw",
    output = "spatLocsObj"
)

del_net <- getSpatialNetwork(
    gobject = viz,
    spat_unit = "aggregate",
    output = "spatialNetworkObj",
    name = "Delaunay_network"
)

knn_net <- getSpatialNetwork(
    gobject = viz,
    spat_unit = "aggregate",
    output = "spatialNetworkObj",
    name = "kNN_network"
)

expr_vals_raw <- getExpression(
    gobject = viz,
    spat_unit = "aggregate",
    values = "raw",
    output = "exprObj"
)

expr_vals_norm <- getExpression(
    gobject = viz,
    spat_unit = "aggregate",
    values = "normalized",
    output = "exprObj"
)

mg_enr <- getSpatialEnrichment(
    gobject = viz,
    name = "cluster_metagene",
    spat_unit = "aggregate",
    output = "spatEnrObj"
)

gpoly <- getPolygonInfo(
    gobject = viz,
    polygon_name = "aggregate",
    return_giottoPolygon = TRUE
)

gpoints <- getFeatureInfo(
    gobject = viz,
    feat_type = "rna",
    return_giottoPoints = TRUE
)

cm <- getCellMetadata(
    gobject = viz,
    spat_unit = "aggregate",
    feat_type = "rna",
    output = "cellMetaObj"
)

fm <- getFeatureMetadata(
    gobject = viz,
    spat_unit = "aggregate",
    feat_type = "rna",
    output = "featMetaObj"
)

dim_red_umap <- getDimReduction(
    gobject = viz,
    spat_unit = "aggregate",
    feat_type = "rna",
    reduction_method = "umap",
    name = "umap",
    output = "dimObj"
)

dim_red_pca <- getDimReduction(
    gobject = viz,
    spat_unit = "aggregate",
    feat_type = "rna",
    reduction_method = "pca",
    name = "pca",
    output = "dimObj"
)

nn <- getNearestNetwork(
    gobject = viz,
    spat_unit = "aggregate",
    feat_type = "rna",
    nn_type = "sNN",
    name = "sNN.pca",
    output = "nnNetObj"
)

viz <- createSpatialGrid(
  gobject = viz,
  sdimx_stepsize = 100,
  sdimy_stepsize = 100
)

spatial_grid_obj <- getSpatialGrid(viz)

# images
# Freshly generate to allow reconnection
data_path <- system.file("/Mini_datasets/Vizgen/Raw/", package = "GiottoData")
DAPI_z0_image_path <- paste0(data_path, "/", "images/mini_dataset_dapi_z0.jpg")
polyT_z0_image_path <- paste0(data_path, "/", "images/mini_dataset_polyT_z0.jpg")

# x and y information from original script
ultra_mini_extent <- terra::ext(c(6400.029, 6900.037, -5150.007, -4699.967))

image_paths <- c(DAPI_z0_image_path, polyT_z0_image_path)
image_names <- c("dapi_z0", "polyT_z0")

imagelist <- createGiottoLargeImageList(
    raster_objects = image_paths,
    names = image_names,
    negative_y = TRUE,
    extent = ultra_mini_extent
)



# 2. save subobjects ####
# ----------------- #

# you need to use your local GiottoData repo
giottodata_repo <- "/Users/gsi-local/Documents/GitHub/GiottoData/"
miniobj_path <- paste0(giottodata_repo, "inst/Mini_objects/subobjects/")

# expression
saveRDS(expr_vals_raw, file = paste0(miniobj_path, "/", "exprObj/viz_agg_expr_raw.RDS"))
saveRDS(expr_vals_norm, file = paste0(miniobj_path, "/", "exprObj/viz_agg_expr_norm.RDS"))

# cell metadata
saveRDS(cm, file = paste0(miniobj_path, "/", "cellMetaObj/viz_agg_cellmeta.RDS"))

# feat metadata
saveRDS(fm, file = paste0(miniobj_path, "/", "featMetaObj/viz_agg_featmeta.RDS"))

# spatial locations
saveRDS(spat_locs, file = paste0(miniobj_path, "/", "spatLocsObj/viz_agg_spatlocs.RDS"))

# spatial network
saveRDS(del_net, file = paste0(miniobj_path, "/", "spatialNetworkObj/viz_agg_spatnet_del.RDS"))
saveRDS(knn_net, file = paste0(miniobj_path, "/", "spatialNetworkObj/viz_agg_spatnet_knn.RDS"))

# spatial enrichment
saveRDS(mg_enr, file = paste0(miniobj_path, "/", "spatEnrObj/viz_agg_metagene.RDS"))

# nearest neighbor network
saveRDS(nn, file = paste0(miniobj_path, "/", "nnNetObj/viz_agg_sNN.RDS"))

# dim reduction
saveRDS(dim_red_pca, file = paste0(miniobj_path, "/", "dimObj/viz_agg_pca.RDS"))
saveRDS(dim_red_umap, file = paste0(miniobj_path, "/", "dimObj/viz_agg_umap.RDS"))

gpoints <- GiottoClass::wrap(gpoints)
saveRDS(gpoints, file = paste0(miniobj_path, "/", "giottoPoints/viz_agg_gpoints.RDS"))

gpoly <- GiottoClass::wrap(gpoly)
saveRDS(gpoly, file = paste0(miniobj_path, "/", "giottoPolygon/viz_agg_gpoly.RDS"))

# largeImages
# loaded by replacing a portion of the path with the end user's library path location
saveRDS(imagelist[[1]], file = paste0(miniobj_path, "/", "giottoLargeImage/viz_dapi_z0.RDS"))
saveRDS(imagelist[[2]], file = paste0(miniobj_path, "/", "giottoLargeImage/viz_polyT_z0.RDS"))

# spatial grid
saveRDS(spatial_grid_obj, file = paste0(miniobj_path, "/", "spatialGridObj/viz_agg_spatialGridObj.RDS"))
