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
store1_simp=gBuffer(store1,width=1,byid=T) # Needed to make geometry valid

#@ Extract projection information from polygon
spatial_proj=proj4string(store1)


# Convert storage polygon to a list of SpatialPolygons (each being a storage area)
store1_list=list()
for(i in 1:length(store1_simp)){
    store1_list[[i]]=SpatialPolygons(store1_simp@polygons[i], proj4string=CRS(spatial_proj))
}


#@ Compute the channel polygon from hec-ras file, and convert it to a spatial format
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))
chan2=gBuffer(chan,width=1, byid=T) # Gets rid of invalid geometry
chan2=SpatialPolygonsDataFrame(chan2, data=data.frame(DN=seq(1,length(chan2))), match.ID=FALSE) # Needed to write to shapefile

#@ Write channel polygon to shapefile for viewing in GIS
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
    for(j in i:length(store1_list)){
        # Note that we only loop from i, so we avoid double counting connections
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


#@
#@ Step 2.1 - Make text describing the storage areas, that can be inserted into hecras .g01 file
#@

storage_text=list()
storage_names=list()
#@ We wish to keep the names of storage areas unique. One way to attempt to
#@ do this is to append a timestamp to the name
name_timestamp=floor(as.numeric(Sys.time())*100)%%1000000+i # Generic stamp for storage name
for(i in 1:length(store1_list)){
    #name=paste('Fake_store',i,sep="") # Name for storage area
    storage_names[[i]]=paste('St',(name_timestamp+i),sep="") # Name for storage area
    storage_text[[i]]=make_storage_area_text(store1_list[[i]], store1_stage_vol_list[[i]], storage_names[[i]])
}

#@ Append to a new hecras file = hec_lines2
hec_lines2=hec_lines
end_storage=grep('Connection=', hec_lines2)[1]-2 # We need to insert storage areas here

hec_tmp=hec_lines2[1:end_storage]
#@ Append storage
for(i in 1:length(storage_text)){
    hec_tmp=c(hec_tmp, " ") # Add a newline
    hec_tmp=c(hec_tmp,storage_text[[i]]) # Add the storage area
}
#@ Append the rest of geometry file
hec_tmp=c(hec_tmp, hec_lines2[(end_storage+1):length(hec_lines2)])
#@ Update the new hecras file
hec_lines2=hec_tmp


#@
#@ Step 2.2 -- Loop over all overlapping storage areas, and make a storage area connection
#@
storage_connection_text_all=c()
for(i in 1:length(storage_intersections)){
  intersections=storage_intersections[[i]]
  if(is.na(intersections)) next
  
  for(j in 1:length(intersections)){
      k=intersections[j]
      storage_connection_text=
          make_storage_connection_text( store1_list[[i]], store1_list[[k]], 
                                        storage_names[[i]], storage_names[[k]],
                                        lidar_DEM )    
      storage_connection_text_all=c(storage_connection_text_all, storage_connection_text, " ")
  }


}

#@
#@ Insert storage area connection text into the output file
#@

#@ find appropriate index
end_storage_inds=max(grep("^Conn HTab HWMax=", hec_lines2)) # Last or second last line of storage areas
if(hec_lines2[end_storage_inds+1]!=""){
    end_storage_inds=end_storage_inds+2
}else{
    end_storage_inds=end_storage_inds+1

}

#@ Do insertion
ll=length(hec_lines2)
hec_linestmp=c(hec_lines2[1:end_storage_inds], 
               storage_connection_text_all, 
               hec_lines2[(end_storage_inds+1):ll])

hec_lines2=hec_linestmp

#@ Write to output
cat(hec_lines2,file='hectest.g05',sep="\n") 
