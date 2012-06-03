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
potential_storage_file='fine_watersheds/watershed.shp'
storage_file_layername='watershed'
dem_filename='mydem.tif'

#@ Read storage polygon
#store1=readOGR(potential_storage_file, layer=storage_file_layername)

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
#chan2=gBoundary(chan2, byid=T)
#writeOGR(chan,dsn='chan_shapefile',layer='chan_shapefile',driver='ESRI Shapefile')

store1_simp=gBuffer(store1,width=1,byid=T) # Needed to make geometry valid
#store1_simp=gBoundary(store1_simp,byid=T)
#store1_simp=gSimplify(store1_simp,tol=30.0, topologyPreserve=T) # Reduce number of polygons
# Try to clip store1_simp to a bounding box
#clipbox=bbox(chan2)
#clipbox_coords=rbind(c(clipbox[1,1], clipbox[2,1]), 
#                     c(clipbox[1,1], clipbox[2,2]),
#                     c(clipbox[1,2], clipbox[2,2]),
#                     c(clipbox[1,2], clipbox[2,1]),
#                     c(clipbox[1,1], clipbox[2,1]))
#clipbox_tmp=SpatialPolygons(list(Polygons(list(Polygon(clipbox_coords)),ID=1)), proj4string=CRS(spatial_proj))               

#store1_simp2=gIntersection(store1_simp,ff_tmp)

#@ Take difference
#store1_nochan=gDifference(store1_simp, chan2,byid=TRUE)


