## R code
library(rgdal)
library(sp)
# Read in utility functions
source('hec_help_util.R')

# Input parameters
hecras_channels_file='May_june_2012.g05'
potential_storage_file='fine_watersheds/watershed.shp'
storage_file_layername='watershed'

# EPSG code of projection system
spatial_epsg=3123
# Convert to proj4string
epsg_codes=make_EPSG()
epsg_relevant=which(epsg_codes[,1]==spatial_epsg)
spatial_proj=epsg_codes[epsg_relevant,3]


# Get channel polygon
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))

# Read storage polygon
store1=readOGR(potential_storage_file, layer=storage_file_layername)


