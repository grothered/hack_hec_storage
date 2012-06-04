#######################################################################
#
#@ R code for semi-automating the creation of storage areas to hec-ras.
#
#######################################################################

library(rgdal)
library(sp)
library(rgeos)

#@ Read in utility functions
source('hec_help_util.R')

#@ Input parameters
hecras_channels_file='May_june_2012.g05'
potential_storage_file='manual_store/storage1.shp'
storage_file_layername='storage1'
dem_filename='mydem.tif'

#@ Read storage polygon
store1=readOGR(potential_storage_file, layer=storage_file_layername)
store1_simp=gBuffer(store1,width=1,byid=T) # Needed to make geometry valid
# Convert to a list of SpatialPolygons
store1_list=list()
for(i in 1:length(store1_simp)){
    store1_list[[i]]=store1_simp@polygons[i]
}

#@ EPSG code of projection system
#spatial_epsg=3123
#@ Convert to proj4string
#epsg_codes=make_EPSG()
#epsg_relevant=which(epsg_codes[,1]==spatial_epsg)
#spatial_proj=epsg_codes[epsg_relevant,3]
spatial_proj=proj4string(store1)

#@ Get channel polygon
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))
chan2=gBuffer(chan,width=1) # Gets rid of invalid geometry
chan2=gSimplify(chan2,tol=1, topologyPreserve=T) # Gets rid of invalid geometries
chan2=SpatialPolygonsDataFrame(chan2, data=data.frame(DN=1), match.ID=FALSE) # Needed to write to shapefile

#@ Write channel polygon to shapefile
writeOGR(chan2,dsn='chan_shapefile3',layer='chan_shapefile',driver='ESRI Shapefile',overwrite=T)

# Convert chan2 to a list of SpatialPolygons
chan2_list=list()
for(i in 1:length(chan2)){
    chan2_list[[i]]=chan2@polygons[i]
}

## Step 1: Compute Stage-Volume relation for every element of store1_list

