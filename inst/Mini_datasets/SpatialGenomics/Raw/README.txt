## Spatial Genomics Mini Kidney Object ####

Raw data for the SG Mini Kidney Object is too large for this repo and can be found at:
https://github.com/drieslab/spatial-datasets/data/2023_spatial_genomics_mouse_kidney/Raw.zip

Use GiottoData's function: getSpatialDataset('sg_mini_kidney') to download the Raw.zip datafile to your working directory.

This folder contains three files: DAPI (cell nuclei coordinates), mask (cell segmentations), and transcripts (coordinates of transcripts and cell number).

The directory to which you save this zip file to can be used to create a Giotto object using the provided script. The datadir variable will automatically set to your working directory OR you can specify a different directory in which your Spatial Genomics files are stored.


##-------------------------------------####
