#######################################################################
#
#@ R code for semi-automating the creation of storage areas in hec-ras.
#
#######################################################################
#
# This code takes as input:
#
# 1) A hec-ras .g0x file (describing the geometry)
#
# 2) A polygon shapefile (storage shapefile) containing several individual
#    polygons which the user would like to add as storage areas to the hec-ras
#    geometry
#
#    If any 2 polygons overlap with each other, it is assumed that a storage area
#    connection is desired at that location
#
#    If any polygon overlaps with a channel in the hec-ras file, it is
#    assumed that a lateral weir connecting the 2 is desired
#
# 3) A raster DEM file
#
#
# It is assumed that 1,2,3 all use the same coordinate system (units in metres)
# , and that 2) and 3) have the correct projection information, and that 1 & 2
# are inside 3.
#
# 
# The code does the following:
#
# 1) Creates a polygon shapefile of the channel network (channel shapefile), by
#    joining boundary points of each cross-section within each river reach in the
#    hec-ras file.
#    Creates a polygon shapefile of the existing storage areas
#    Creates a polygon shapefile of the channel bank points
#
#    The code will stop here if 'create_shapefiles_of_existing_rasfile=TRUE', which can be useful if you
#    just want to create a channel polygon
#
# 2) For each individual polygon in the storage shapefile, add it as a storage area in
#    the hec-ras geometry file, using a stage-volume relation computed from the
#    raster inside the polygon. 
#
# 3) Compute all intersections of polygons in the storage-shapefile. These are taken as
#    sites where a 'storage area connection' is to be located. The connection
#    is created as a weir, with the heights over the weir reflecting the
#    heights in the DEM within the intersection zone. This is added to the
#    hec-ras geometry file.
#    [Note that the weir length = intersection_polygon_boundary_length/3, which
#    should be approximately correct for a relatively wide, narrow intersection
#    polygon. This could potentially be done in a more refined way].
#
# 4) Computes all the intersections of polygons in the storage shapefile, with
#    the channel shapefile These are taken as sites where a 'lateral weir
#    connection' is to be located. The elevation of the weir is selected based on
#    the elevations of the raster DEM over that zone.
#
#######################################################################

#@ Input parameters
hecras_channels_file='May_june_2012.g05'
output_file='hectest.g05'
potential_storage_file='manual_store/storage_2.shp'
lidar_DEM_file='C:/Users/Gareth/Documents/work/docs/Nov_2011_workshops/qgis/LIDAR_and_IMAGERY/DEM/10m_DEM/test2_10m.tif'
#lidar_DEM_file='/media/Windows7_OS/Users/Gareth/Documents/work/docs/Nov_2011_workshops/qgis/LIDAR_and_IMAGERY/DEM/10m_DEM/test2_10m.tif'
create_shapefiles_of_existing_rasfile=TRUE
vertical_datum_offset=10.5


storage_file_layername=basename(potential_storage_file) # Might need to edit this
storage_file_layername=strsplit(storage_file_layername, "\\.")[[1]][1]

#@ Get libraries
library(rgdal)
library(sp)
library(rgeos)
library(raster)

#@ Read in utility functions
source('hec_help_util.R')

# Read lidar_dem
lidar_DEM=raster(lidar_DEM_file)

#@ Extract projection information from lidar
spatial_proj=lidar_DEM@crs@projargs

#@ Compute the channel polygon from hec-ras file, and convert it to a spatial format
print('EXTRACTING CHANNEL POLYGON')
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))
chan2=gBuffer(chan,width=0.1, byid=T) # Gets rid of invalid geometry
chan2=SpatialPolygonsDataFrame(chan2, data=data.frame(DN=seq(1,length(chan2))), match.ID=FALSE) # Needed to write to shapefile

print("EXTRACTING CHANNEL XSECT BOUNDARY POINTS")
chan_boundary_points=make_channel_boundary_points(hecras_channels_file, spatial_proj)

print('EXTRACTING EXISTING STORAGE AREAS FROM HECRAS FILE')
print(' ')
old_storage=get_existing_storage_areas(hecras_channels_file,spatial_proj)

if(create_shapefiles_of_existing_rasfile){
    #@ Write channel polygon to shapefile for viewing in GIS
    try(writeOGR(chan2,dsn='chan_shapefile',layer='chan_shapefile',driver='ESRI Shapefile',overwrite=TRUE))

    #@ Write channel boundary points to a shapefile for viewing in GIS
    try(writeOGR(chan_boundary_points,dsn='chan_points',layer='chan_points',driver='ESRI Shapefile',overwrite=TRUE))

    #@ Write existing storage areas to shapefile
    try(writeOGR(old_storage, dsn='existing_storage',layer='existing_storage',driver='ESRI Shapefile',overwrite=TRUE))
    
    print(' ')
    print('Tried to create polygons from hec-ras files creation') 
    print('I will not proceed further while create_shapefiles_of_existing_rasfile=TRUE.')
    print('If you got a "Layer creation failed" error, it probably means')
    print(' that one of the shapefiles is already open,')
    print(' in which case R may not be allowed to update it')

}else{

    #@ Read storage polygon
    new_store=readOGR(potential_storage_file, layer=storage_file_layername)
    new_store_simp=gBuffer(new_store,width=0.1,byid=T) # Helpful to make geometry valid


    ############################ FORMAT CONVERSION for R

    #@ Convert new storage polygon to a list of SpatialPolygons (each being a storage area)
    new_store_list=list()
    for(i in 1:length(new_store_simp)){
        new_store_list[[i]]=SpatialPolygons(new_store_simp@polygons[i], proj4string=CRS(spatial_proj))
    }

    old_store_list=list()
    for(i in 1:length(old_storage@polygons)){
        old_store_list[[i]]=SpatialPolygons(old_storage@polygons[i], proj4string=CRS(spatial_proj))
    }

    #@ Convert chan2 to a list of SpatialPolygons
    chan2_list=list()
    for(i in 1:length(chan2)){
        chan2_list[[i]]=SpatialPolygons(chan2@polygons[i], proj4string=CRS(spatial_proj))
    }

    ########################### MAIN CODE

    #@ Step 1: Compute Stage-Volume relation for every element of new_store_list
    print("COMPUTING STAGE-VOLUME RELATIONS FOR STORAGE AREAS")
    new_store_stage_vol_list=list()
    for(i in 1:length(new_store_list)){
        my_poly=new_store_list[[i]]
        new_store_stage_vol_list[[i]] = compute_stage_vol_relation(my_poly,lidar_DEM, vertical_datum_offset)
    }

    #@ Step 1.5: Identify intersections of storage polygons with each other, or with the channel network
    print("COMPUTING 'NEW STORAGE AREA' INTERSECTIONS")
    new_storage_intersections=list()
    for(i in 1:length(new_store_list)){
        intersections=c()
        for(j in i:length(new_store_list)){
            # Note that we only loop from i, so we avoid double counting connections
            if(i==j) next

            if(gIntersects(new_store_list[[i]], new_store_list[[j]])){
                intersections=c(intersections,j)
            } 
        }
        if(is.null(intersections)){
            new_storage_intersections[[i]]=NA
        }else{
            new_storage_intersections[[i]]=intersections
        }
    } 
    
    print("COMPUTING 'NEW STORAGE' - 'OLD STORAGE' AREA INTERSECTIONS")
    old_storage_intersections=list()
    if(length(old_store_list)>0){
        for(i in 1:length(new_store_list)){
            intersections=c()
            for(j in 1:length(old_store_list)){
                
                if(gIntersects(new_store_list[[i]], old_store_list[[j]])){
                    intersections=c(intersections,j)
                } 
            }
            if(is.null(intersections)){
                old_storage_intersections[[i]]=NA
            }else{
                old_storage_intersections[[i]]=intersections
            }
        }
    }


    print("COMPUTING 'NEW STORAGE' - CHANNEL AREA INTERSECTIONS")
    channel_intersections=list()
    for(i in 1:length(new_store_list)){
        intersections=c()
        for(j in 1:length(chan2_list)){
            if(gIntersects(new_store_list[[i]], chan2_list[[j]])){
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
    print("CONSTRUCTING 'NEW STORAGE AREA' GEOMETRIES")
    storage_text=list()
    new_storage_names=list()
    #@ We wish to keep the names of storage areas unique. One way to attempt to
    #@ do this is to append a timestamp to the name
    name_timestamp=floor(as.numeric(Sys.time())*100)%%1000000+i # Generic stamp for storage name
    for(i in 1:length(new_store_list)){
        #name=paste('Fake_store',i,sep="") # Name for storage area
        new_storage_names[[i]]=paste('St',(name_timestamp+i),sep="") # Name for storage area
        storage_text[[i]]=make_storage_area_text(new_store_list[[i]], new_store_stage_vol_list[[i]], new_storage_names[[i]])
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
    print("CONSTRUCTING STORAGE AREA CONNECTION GEOMETRIES AT 'NEW STORAGE'-'NEW STORAGE' INTERSECTIONS")
    storage_connection_text_all=c()
    for(i in 1:length(new_storage_intersections)){
      intersections=new_storage_intersections[[i]]
      if(is.na(intersections)) next
      
      for(j in 1:length(intersections)){
          k=intersections[j]
          storage_connection_text=
              make_storage_connection_text( new_store_list[[i]], new_store_list[[k]], 
                                            new_storage_names[[i]], new_storage_names[[k]],
                                            lidar_DEM,vertical_datum_offset)    
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


    #@
    #@ Step 2.3 -- Loop over all overlapping storage areas, and make a storage area connection
    #@
    print("CONSTRUCTING STORAGE AREA CONNECTION GEOMETRIES AT 'NEW STORAGE'-'OLD STORAGE' INTERSECTIONS")
    storage_connection_text_all=c()
    for(i in 1:length(old_storage_intersections)){
      intersections=old_storage_intersections[[i]]
      if(is.na(intersections)) next
      
      for(j in 1:length(intersections)){
          k=intersections[j]
          storage_connection_text=
              make_storage_connection_text( new_store_list[[i]], old_store_list[[k]], 
                                            new_storage_names[[i]], old_storage$names[k],
                                            lidar_DEM,vertical_datum_offset)    
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



    #@
    #@
    #@ Make a spatial points object holding the channel cross-section end points, and their section label
    #@ We can use this to associate lateral structures with x-sections 
    #@
    #@
    #stop()
    # Loop over all storage-channel intersections
    print("CREATING LATERAL STRUCTURES AT 'NEW STORAGE'-CHANNEL INTERSECTIONS")
    hec_linestmp=hec_lines2 # Copy output file for modification
    for(i in 1:length(channel_intersections)){

        if(is.na(channel_intersections[[i]])) next

        for(j in 1:length(channel_intersections[[i]])){
            k=channel_intersections[[i]][j]
            #@ Iteratively update hec_linestmp by inserting laterl weir
            hec_linestmp=make_lateral_weir_text(new_store_list[[i]], new_storage_names[[i]],
                                                     chan2_list[[k]], chan_boundary_points, 
                                                     lidar_DEM, vertical_datum_offset, 
                                                     hec_linestmp)

        }

    }

    hec_lines2=hec_linestmp
    #@ Write to output
    cat(hec_lines2,file=output_file,sep="\n") 

}
