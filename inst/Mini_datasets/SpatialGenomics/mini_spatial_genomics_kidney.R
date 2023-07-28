## Mini Spatial Genomics Kidney script and dataset preparation##
# ------------------------- #
library(Giotto)

# 1. Object Creation & Filtering #
# ------------------------- #
# Set directory containing SG data
datadir = 'secret/path/to/Amelia/directories'

# Create SG object using function
sg <- createSpatialGenomicsObject(sg_dir = datadir)

# Aggregate
sg = calculateOverlapRaster(sg, spatial_info = 'cell', feat_info = 'rna')
sg = overlapToMatrix(sg)
sg = addSpatialCentroidLocations(sg)

# Normalize
sg = filterGiotto(sg, feat_det_in_min_cells = 10, min_det_feats_per_cell = 2, expression_threshold = 1)
sg = normalizeGiotto(sg)

# Add statistics
sg = addStatistics(sg)

# Total transcripts after filtering
# low, mid, high
custom_scale = c('#440154', '#1F968B', '#FDE725')

spatInSituPlotPoints(sg,
                     show_polygon = TRUE,
                     polygon_alpha = 1,
                     polygon_fill_gradient = custom_scale,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE,
                     polygon_line_size = 0.1)

# 2. Dimension Reduction #
# ------------------------- #
# Calculate Highly Variable Features
sg = calculateHVF(gobject = sg)
cat(fDataDT(sg)[, sum(hvf == 'yes')], 'hvf found')
# hvf = 18 -> better to use all genes -> feats_to_use = NULL
sg = runPCA(gobject = sg,
            spat_unit = 'cell',
            expression_values = 'scaled',
            feats_to_use = NULL,
            scale_unit = F,
            center = F)

# Visualize Screeplot and PCA
screePlot(sg,
          ncp = 20,
          save_param = list(
            save_name = 'sg_screePlot'))
showGiottoDimRed(sg)
plotPCA(sg,
        spat_unit = 'cell',
        dim_reduction_name = 'pca',
        dim1_to_use = 1,
        dim2_to_use = 2)

# Run and Plot tSNE and UMAP
sg = runtSNE(sg,
             dimensions_to_use = 1:10,
             spat_unit = 'cell',
             check_duplicates = FALSE)
sg = runUMAP(sg,
             dimensions_to_use = 1:10,
             spat_unit = 'cell')
plotTSNE(sg,
         point_size = 0.5,
         save_param = list(
           save_name = 'sg_tSNE'))
plotUMAP(sg,
         point_size = 0.5,
         save_param = list(
           save_name = 'sg_UMAP'))

# 3. Clustering #
# ------------------------- #
# Leiden clustering and UMAP cluster visualization
sg = createNearestNetwork(sg,
                          dimensions_to_use = 1:10,
                          k = 10,
                          spat_unit = 'cell')
sg = doLeidenCluster(sg,
                     resolution = 0.25,
                     n_iterations = 100,
                     spat_unit = 'cell')
plotUMAP(gobject = sg,
         spat_unit = 'cell',
         cell_color = 'leiden_clus',
         show_legend = FALSE,
         point_size = 0.5,
         point_shape = 'no_border')
# Plot Leiden clusters onto spatial image plot and as InSitu Polygons
my_spatPlot <- spatPlot2D(gobject = sg,
                          spat_unit = 'cell',
                          cell_color = 'leiden_clus',
                          point_size = 0.8,
                          point_shape = 'no_border',
                          show_legend = TRUE)
spatInSituPlotPoints(sg,
                     show_image = FALSE,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_line_size = 0.2,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE)

# 4. Identify Spatial Genes #
# ------------------------- #
# Establish Delaunay Network
sg = createSpatialNetwork(gobject = sg, minimum_k = 2,
                          maximum_distance_delaunay = 250)
sg = createSpatialNetwork(gobject = sg, minimum_k = 2,
                          method = 'kNN', k = 10)
spatPlot(gobject = sg, show_network = T,
         network_color = 'blue', spatial_network_name = 'Delaunay_network',
         point_size = 1.0, cell_color = 'leiden_clus')

# By k-means
km_spatialgenes = binSpect(sg)
spatFeatPlot2D(sg, expression_values = 'scaled',
               feats = km_spatialgenes[1:4]$feats,
               point_shape = 'border', point_border_stroke = 0.1,
               show_network = F, network_color = 'lightgrey', point_size = 1.0,
               cow_n_col = 2)

# 5. Save Giotto Object #
# ------------------------- #
format(object.size(sg), units = 'Mb')

# Use your local GiottoData repo
giottodata_repo = 'secret/path/to/Amelia/Mini_datasets_directories/'

saveGiotto(gobject = sg,
           foldername = 'SpatialGenomicsObject',
           dir = paste0(giottodata_repo, '/', 'SpatialGenomics/'),
           overwrite = TRUE)

# Test by Loading Object 
sg_object_dir = 'secret/path/to/Amelia/Mini_datasets_directories/SpatialGenomics/SpatialGenomicsObject/'
get_sg = loadGiotto(path_to_folder = sg_object_dir)
spatInSituPlotPoints(get_sg,
                     show_image = FALSE,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_alpha = 1,
                     polygon_line_size = 0.2,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE,
                     coord_fix_ratio = TRUE)
