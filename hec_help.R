#######################################################################
#
#@ R code for semi-automating the creation of storage areas to hec-ras.
#
#######################################################################

library(rgdal)
library(sp)
library(rgeos)
library(raster)

#@ Read in utility functions
source('hec_help_util.R')

#@ Input parameters
hecras_channels_file='May_june_2012.g05'
potential_storage_file='manual_store/storage1.shp'
storage_file_layername='storage1'
lidar_DEM_file='C:/Users/Gareth/Documents/work/docs/Nov_2011_workshops/qgis/LIDAR_and_IMAGERY/DEM/10m_DEM/test2_10m.tif'

# Read lidar_dem
lidar_DEM=raster(lidar_DEM_file)

#@ Read storage polygon
store1=readOGR(potential_storage_file, layer=storage_file_layername)
spatial_proj=proj4string(store1)
store1_simp=gBuffer(store1,width=1,byid=T) # Needed to make geometry valid
# Convert to a list of SpatialPolygons
store1_list=list()
for(i in 1:length(store1_simp)){
    store1_list[[i]]=SpatialPolygons(store1_simp@polygons[i], proj4string=CRS(spatial_proj))
}


#@ Get channel polygon from hec-ras file, and convert to spatial format
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))
chan2=gBuffer(chan,width=1, byid=T) # Gets rid of invalid geometry
#chan2=gSimplify(chan2,tol=1, topologyPreserve=T) # Gets rid of invalid geometries
chan2=SpatialPolygonsDataFrame(chan2, data=data.frame(DN=seq(1,length(chan2))), match.ID=FALSE) # Needed to write to shapefile

#@ Write channel polygon to shapefile
writeOGR(chan2,dsn='chan_shapefile5',layer='chan_shapefile',driver='ESRI Shapefile',overwrite=TRUE)

#@ Convert chan2 to a list of SpatialPolygons
chan2_list=list()
for(i in 1:length(chan2)){
    chan2_list[[i]]=SpatialPolygons(chan2@polygons[i], proj4string=CRS(spatial_proj))
}

#@ Step 1: Compute Stage-Volume relation for every element of store1_list
print("COMPUTING STAGE-VOLUME RELATIONS FOR STORAGE AREAS")
store1_stage_vol_list=list()
for(i in 1:length(store1_list)){
    my_poly=store1_list[[i]]
    store1_stage_vol_list[[i]] = compute_stage_vol_relation(my_poly,lidar_DEM, vertical_datum_offset=10.5)
}

#@ Step 1.5: Identify intersections of storage polygons with each other, or with the channel network
print("COMPUTING STORAGE AREA INTERSECTIONS")
storage_intersections=list()
for(i in 1:length(store1_list)){
    intersections=c()
    for(j in 1:length(store1_list)){
        if(i==j) next

        if(gIntersects(store1_list[[i]], store1_list[[j]])){
            intersections=c(intersections,j)
        } 
    }
    if(is.null(intersections)){
        storage_intersections[[i]]=NA
    }else{
        storage_intersections[[i]]=intersections
    }
} 

print("COMPUTING CHANNEL AREA INTERSECTIONS")
channel_intersections=list()
for(i in 1:length(store1_list)){
    intersections=c()
    for(j in 1:length(chan2_list)){
        if(gIntersects(store1_list[[i]], chan2_list[[j]])){
            intersections=c(intersections,j)
        } 
    }
    if(is.null(intersections)){
        channel_intersections[[i]]=NA
    }else{
        channel_intersections[[i]]=intersections
    }
} 

#@ Step 2: Make a new hec-ras geometry file, which includes the rivers /
#@ junctions and storage areas only?  Then we can add new storage area
#@ connections / lateral weirs to this. Then, the user can import these only,
#@ using hecs import functionality. This seems safer than trying to
#@ automatically avoid duplicating / skipping already defined storage areas

 
#@ Read file
fin=file(hecras_channels_file, open='r')
hec_lines=readLines(fin)
close(fin)

#@ Make text describing the storage areas, that can be inserted into hecras .g01 file
storage_text=list()
for(i in 1:length(store1_list)){
    name=paste('Fake_store',i,sep="")
    storage_text[[i]]=make_storage_area_text(store1_list[[i]], store1_stage_vol_list[[i]], name)
}

#@ Append to a new hecras file
hec_lines2=hec_lines
end_storage=grep('Connection=', hec_lines2)[1]-2 # We need to insert storage areas here

hec_tmp=hec_lines2[1:end_storage]
#@ Append storage
for(i in 1:length(storage_text)){
    hec_tmp=c(hec_tmp, " ")
    hec_tmp=c(hec_tmp,storage_text[[i]])
}
#@ Append rest of file
hec_tmp=c(hec_tmp, hec_lines2[(end_storage+1):length(hec_lines2)])

hec_lines2=hec_tmp

cat(hec_lines2,file='hectest.g05',sep="\n") 
