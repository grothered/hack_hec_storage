#######################################################################
#
#@ R code for semi-automating the creation of storage areas in hec-ras.
#@ , and for doing some other useful things
#
#######################################################################
#
# This code takes as input:
#
# 1) A hec-ras .g0x file (describing the geometry), WHERE ALL THE CROSS-SECTIONAL CUTLINES ARE GEOREFERENCED
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
# The code does the following, depending on how you set it up
#
# 1) Creates a polygon shapefile of the channel network (channel shapefile), by
#    joining boundary points of each cross-section within each river reach in the
#    hec-ras file.
#    Creates a polygon shapefile of the existing storage areas
#    Creates a polygon shapefile of the channel bank points
#
#    These files are written to hard-coded directories in the local directory
#
#    THE CODE WILL STOP HERE if 'create_shapefiles_of_existing_rasfile=TRUE', which can be useful if you
#    just want to create a channel polygon
#
# 1.5) Use the channel network and xsect cutlines to re-define the 'downstream_distance' values in the hec-ras file.
#      Only activated if update_downstream_dist=TRUE. THE CODE WILL STOP HERE IF THIS OPTION IS ON.
#      It computes the channel downstream distance 'along the channel', and the overbank distances are set equal to 
#      the 'straight line distance' between the channel intersection points of consecutive cross-sections
#      WARNING: Make sure you check the distances where the river does funny things [e.g. overlaps a cross-section in multiple places]
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
#
#    [Note that the weir length = intersection_polygon_boundary_length/3, which
#    should be approximately correct for a relatively wide, narrow intersection
#    polygon. This could potentially be done in a more refined way].
#
#    [Note that at present, if the weir length is computed as < 5m, then no weir is constructed
#     This is to prevent very small intersections producing a storage area connection]
#
# 4) Computes all the intersections of polygons in the storage shapefile, with
#    the channel shapefile.
#
#    These are taken as sites where a 'lateral weir connection' is to be
#    located, SO LONG AS THEY CONTAIN MORE THAN 2 BANK POINTS. 
#
#    The elevation of the weir is selected based on the elevations of the
#    raster DEM over that zone.
#
#    NOTE: Because hec-ras cannot have overlapping lateral weirs, we only
#    connect the storage area to bank points which are not in any other storage
#    area. 
#
#    NOTE: Hec-RAS can only have weirs covering < 100 cross-sections. We
#    truncate weirs to meet this requirement
#
#######################################################################

#@ Input parameters
hecras_channels_file='May_june_2012_c.g03' #'south_of_pasig_merg.g01' #'north_of_pasig_merg.g02' #'May_june_2012.g27' 
#hecras_channels_file='/media/Windows7_OS/Users/Gareth/Documents/work/docs/may_june2012_workshops/hec_ras/model/store_geometry_and_boundaries/2012_07_17/May_june_2012.g29'
output_file='hectest.g05'

#potential_storage_file='storage_areas_gt_100x100b/storage_areas_gt_100x100b.shp'
potential_storage_file='C:/Users/Gareth/Documents/work/docs/may_june2012_workshops/hec_ras/hec_helper/new_storeb/newstore.shp'

lidar_DEM_file='C:/Users/Gareth/Documents/work/docs/Nov_2011_workshops/qgis/LIDAR_and_IMAGERY/DEM/10m_DEM/test2_10m.tif'
#lidar_DEM_file='/media/Windows7_OS/Users/Gareth/Documents/work/docs/Nov_2011_workshops/qgis/LIDAR_and_IMAGERY/DEM/10m_DEM/test2_10m.tif'
lidar_DSM_file="F:/manila_DSM/Tiles1km_1km/manila_1m_dsm.vrt"
vertical_datum_offset=10.5


limit_weir_elevation_by_channel_bank_elevation=FALSE # If TRUE, then we lateral weir elevation is forced >= the elevation of the nearby channel bank points. Experimentation suggests it is better to be FALSE.
lower_limit_on_lateral_weir_elevations=-Inf # The weir elevation in hec-ras will always be >= lower_limit_on_lateral_weir_elevation. Set to -Inf to have no limit
logfile='Rlog.log'

#@ Options to do different things
create_shapefiles_of_existing_rasfile=TRUE # IF TRUE, then make shapefiles of the channel geometry, and stop
update_downstream_dist=FALSE # If true, then use the cutlines and channel centrelines in the hecras file to update the downstream distance parameters, then stop 


###################################################################################################


# Start sending output to a file
sink(logfile, split=TRUE)

# Separate filename and basename in the storage_area_file
storage_file_layername=basename(potential_storage_file) # Might need to edit this
storage_file_layername=strsplit(storage_file_layername, "\\.")[[1]][1]

#@ Get libraries
library(rgdal)
library(sp)
library(rgeos)
library(raster)

#@ Read in utility functions
source('hec_help_util.R')

#@ Read lidar_dem
lidar_DEM=raster(lidar_DEM_file)
lidar_DSM=raster(lidar_DSM_file)

#@ Extract projection information from lidar
spatial_proj=lidar_DEM@crs@projargs

#@ Compute the channel polygon from hec-ras file, and convert it to a spatial format
print('EXTRACTING CHANNEL POLYGON')
chan=create_channel_polygon(hecras_channels_file, CRS(spatial_proj))
chan2=gBuffer(chan,width=0.1, byid=T) # Gets rid of invalid geometry
chan2=SpatialPolygonsDataFrame(chan2, data=data.frame(DN=seq(1,length(chan2))), match.ID=FALSE) # Needed to write to shapefile

print("EXTRACTING CHANNEL XSECT BOUNDARY POINTS")
chan_boundary_points=make_channel_boundary_points(hecras_channels_file, spatial_proj)

print("EXTRACTING CHANNEL CUTLINES")
chan_cutlines=make_channel_cutlines(hecras_channels_file, spatial_proj)

print('COMPUTING CENTRELINES')
centrelines=compute_centrelines(hecras_channels_file, spatial_proj)

print('EXTRACTING EXISTING STORAGE AREAS FROM HECRAS FILE')
print(' ')
old_storage=get_existing_storage_areas(hecras_channels_file,spatial_proj)

if(create_shapefiles_of_existing_rasfile){
    #@ Write channel polygon to shapefile for viewing in GIS
    try(writeOGR(chan2,dsn='chan_shapefile',layer='chan_shapefile',driver='ESRI Shapefile',overwrite=TRUE))

    #@ Write channel boundary points to a shapefile for viewing in GIS
    try(writeOGR(chan_boundary_points,dsn='chan_points',layer='chan_points',driver='ESRI Shapefile',overwrite=TRUE))

    #@ Write channel cutlines to a shapefile
    try(writeOGR(chan_cutlines,dsn='chan_cutlines',layer='chan_cutlines',driver='ESRI Shapefile',overwrite=TRUE))
    
    #@ Write channel centrelines to a shapefile
    try(writeOGR(centrelines,dsn='chan_centrelines',layer='centrelines',driver='ESRI Shapefile',overwrite=TRUE))

    #@ Write existing storage areas to shapefile
    if(!is.null(old_storage)){
        try(writeOGR(old_storage, dsn='existing_storage',layer='existing_storage',driver='ESRI Shapefile',overwrite=TRUE))
    } 
    print(' ')
    print('Hopefully I created shapefiles depicting the hec-ras channels') 
    print('I will not proceed further while create_shapefiles_of_existing_rasfile=TRUE.')
    print('If you got a \"Layer creation failed\" error, it probably means')
    print(' that one of the shapefiles is already open, and you are using windows,')
    print(' in which case R may not be allowed to update it')
    print(' In that case, close the shapefile, and try again')

    sink() # Close the file sink
}else{

    #@ Read hecras file
    fin=file(hecras_channels_file, open='r')
    hec_lines=readLines(fin)
    close(fin)

    #@ Get Line numbers associated with the bridges -- we need this much later
    bridge_lines=grep('Bridge Culvert', hec_lines) 

    #@ Adjust the downstream distances based on the spatial data in the hecras file
    if(update_downstream_dist){
        print("COMPUTING UPDATED DOWNSTREAM DISTANCES")
        hec_lines = update_downstream_distances(hec_lines, chan_cutlines, centrelines)
        #browser()
        cat(hec_lines,file=output_file,sep="\n") 
        sink()
        stop('Will not continue until update_downstream_dist=FALSE')
    }

    #@ Read storage polygon
    new_store=readOGR(potential_storage_file, layer=storage_file_layername)
    new_store_simp=gBuffer(new_store,width=0.3,byid=T) # Helpful to make geometry valid, and to make 'just touching' polygons intersect


    ############################ FORMAT CONVERSION for R

    #@ Convert new storage polygon to a list of SpatialPolygons (each being a storage area)
    new_store_list=list()
    for(i in 1:length(new_store_simp)){
        new_store_list[[i]]=SpatialPolygons(new_store_simp@polygons[i], proj4string=CRS(spatial_proj))
    }

    # FIXME: HACK For manila hec-ras work in the following if statement
    if(TRUE){
        print('Removing overly large intersections from the storage area files')
        # Remove intersections that are too large from new_store_list, by merging polygons
        # Note that intersections will occur at manually defined boundaries,
        # so doing this will still respect catchment type boundaries
        new_store_list_less_intersect=list()
        counter=0
        for(i in 1:(length(new_store_list)-1)){
            APPEND=TRUE
            for(j in (i+1):length(new_store_list)){

                if(gIntersects(new_store_list[[i]], new_store_list[[j]])){
                    # Check how large the intersection is compared with each polygon
                    ij_intersect=gIntersection(new_store_list[[i]], new_store_list[[j]])
                    if(gArea(ij_intersect) > 0.2*min(gArea(new_store_list[[i]]), gArea(new_store_list[[j]]))){
                        # Intersection is too large -- add polygon i to j, and move on
                        #browser()
                        new_store_list[[j]]=gUnion(new_store_list[[i]], new_store_list[[j]])
                        APPEND=FALSE
                    }
                }
                
            }
            # If there is no overly large intersection, keep
            if(APPEND==TRUE){
                counter=counter+1
                new_store_list_less_intersect[[counter]]=new_store_list[[i]]
            }
        }
        # Add in the last polygon in new_store_list    
        counter=counter+1
        new_store_list_less_intersect[[counter]] = new_store_list[[length(new_store_list)]]

        # Update new_store_list
        new_store_list=new_store_list_less_intersect
    }# End HACK FOR MANILA WORK
    #browser()

    # Buffer
    for(i in 1:length(new_store_list)){
        new_store_list[[i]] = gBuffer(new_store_list[[i]], width=1.0)
    }

    if(!is.null(old_storage)){
        old_store_list=list()
        for(i in 1:length(old_storage@polygons)){
            old_store_list[[i]]=SpatialPolygons(old_storage@polygons[i], proj4string=CRS(spatial_proj))
        }
    }

    #@ Convert chan2 to a list of SpatialPolygons
    chan2_list=list()
    for(i in 1:length(chan2)){
        chan2_list[[i]]=SpatialPolygons(chan2@polygons[i], proj4string=CRS(spatial_proj))
    }

    ########################### MAIN CODE

    #@ Step 1.: Identify intersections of storage polygons with each other, or with the channel network
    #@           While we do this, we make a set of polygons including the unique parts of storage areas
    print("COMPUTING 'NEW STORAGE'-'NEW STORAGE' AREA INTERSECTIONS")
    new_storage_intersections=list()
    unique_new_store_list=list() # This will store the unique part of the storage area
    for(i in 1:length(new_store_list)){
        #unique_new_store_list[[i]]=gSimplify(new_store_list[[i]], tol=0.1, topologyPreserve=TRUE) # Trick to try to force valid geometries
        unique_new_store_list[[i]]=new_store_list[[i]]

        intersections=c()
        for(j in i:length(new_store_list)){
            # Note that we only loop from i, so we avoid double counting connections
            if(i==j) next

            if(gIntersects(unique_new_store_list[[i]], new_store_list[[j]])){
                intersections=c(intersections,j)
                if(FALSE){
                    # Record the non-intersecting area -- can be prone to problems
                    unique_new_store_list[[i]] = gDifference(unique_new_store_list[[i]], new_store_list[[j]])
                    unique_new_store_list[[i]] = gBuffer(unique_new_store_list[[i]], width=0.) # Hack to make valid
                }
            } 
        }
        if(is.null(intersections)){
            new_storage_intersections[[i]]=NA
        }else{
            new_storage_intersections[[i]]=intersections
        }
    } 
   
    if(!is.null(old_storage)){ 
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
    
    #@ Step 1.5: Compute Stage-Volume relation for every element of new_store_list
    print("COMPUTING STAGE-VOLUME RELATIONS FOR STORAGE AREAS")
    new_store_stage_vol_list=list()
    for(i in 1:length(new_store_list)){
        my_poly=new_store_list[[i]]
        new_store_stage_vol_list[[i]] = compute_stage_vol_relation(my_poly,lidar_DEM, vertical_datum_offset)
    }


    #@ Step 2: Make a new hec-ras geometry file, which includes the rivers /
    #@ junctions and storage areas only?  Then we can add new storage area
    #@ connections / lateral weirs to this. Then, the user can import these only,
    #@ using hecs import functionality. 
     

    #@
    #@ Step 2.1 - Make text describing the storage areas, that can be inserted into hecras .g01 file
    #@
    print("CONSTRUCTING 'NEW STORAGE AREA' GEOMETRIES")
    storage_text=list()
    new_storage_names=list()
    #@ We wish to keep the names of storage areas unique. One way to attempt to
    #@ do this is to append a timestamp to the name
    name_timestamp=floor(as.numeric(Sys.time())*100)%%10000000+i # Generic stamp for storage name
    for(i in 1:length(new_store_list)){
        #name=paste('Fake_store',i,sep="") # Name for storage area
        new_storage_names[[i]]=paste('S',(name_timestamp+i),sep="") # Name for storage area
        storage_text[[i]]=make_storage_area_text(new_store_list[[i]], new_store_stage_vol_list[[i]], new_storage_names[[i]])
    }

    #@ Append to a new hecras file = hec_lines2
    hec_lines2=hec_lines
    end_storage=grep('Connection=', hec_lines2)[1]-2 # We need to insert storage areas here
    if(is.na(end_storage)){
        end_storage=grep('Chan Stop Cuts', hec_lines2)[1] -2 
    }

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

    #save(hec_lines2, file='tmp_new_storage_geo.Rdata')
    save.image(file='workspace_keep_new_storage_geo.Rdata')
    #@
    #@ Step 2.2 -- Loop over all overlapping storage areas, and make a storage area connection
    #@
    print("CONSTRUCTING STORAGE AREA CONNECTION GEOMETRIES AT 'NEW STORAGE'-'NEW STORAGE' INTERSECTIONS")
    storage_connection_text_all=c()
    for(i in 1:length(new_storage_intersections)){
      intersections=new_storage_intersections[[i]]
      if(is.na(intersections[1])) next
      
      for(j in 1:length(intersections)){
          k=intersections[j]
          storage_connection_text=
              make_storage_connection_text( new_store_list[[i]], new_store_list[[k]], 
                                            new_storage_names[[i]], new_storage_names[[k]],
                                            lidar_DEM,vertical_datum_offset)    
          if(!is.na(storage_connection_text[1])){
              storage_connection_text_all=c(storage_connection_text_all, storage_connection_text, " ")
          }
      }


    }

    #@
    #@ Insert storage area connection text into the output file
    #@

    #@ find appropriate index
    end_storage_inds=max(grep("^Conn HTab HWMax=", hec_lines2)) # Last or second last line of storage areas
    if(!is.finite(end_storage_inds)){
        end_storage_inds=grep('Chan Stop Cuts', hec_lines2)[1] -2 
    }

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

    #save(hec_lines2, file='tmp_storage_con_A.Rdata')
    save.image('workspace_keep_storage_con_A.Rdata')

    #@
    #@ Step 2.3 -- Loop over all overlapping storage areas, and make a storage area connection
    #@
    if(!is.null(old_storage)){
        print("CONSTRUCTING STORAGE AREA CONNECTION GEOMETRIES AT 'NEW STORAGE'-'OLD STORAGE' INTERSECTIONS")
        storage_connection_text_all=c()
        for(i in 1:length(old_storage_intersections)){
          intersections=old_storage_intersections[[i]]
          if(is.na(intersections[1])) next
          
          for(j in 1:length(intersections)){
              k=intersections[j]
              storage_connection_text=
                  make_storage_connection_text( new_store_list[[i]], old_store_list[[k]], 
                                                new_storage_names[[i]], old_storage$name[k],
                                                lidar_DEM,vertical_datum_offset)    
              if(!is.na(storage_connection_text[1])){
                  storage_connection_text_all=c(storage_connection_text_all, storage_connection_text, " ")
              }
          }
        }
    }

    #@
    #@ Insert storage area connection text into the output file
    #@

    #@ find appropriate index
    end_storage_inds=max(grep("^Conn HTab HWMax=", hec_lines2)) # Last or second last line of storage areas
    if(!is.finite(end_storage_inds)){
        end_storage_inds=grep('Chan Stop Cuts', hec_lines2)[1] -2 
    }
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


    #save(hec_lines2, file='tmp_storage_con_B.Rdata')
    save.image('workspace_keep_storage_con_B.Rdata')

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

        if(is.na(channel_intersections[[i]][1])) next
            
        #@ First:
        unique_new_store_list[[i]]=gBuffer(unique_new_store_list[[i]],width=2.) # Buffer helps to avoid skipping the odd point

        for(j in 1:length(channel_intersections[[i]])){
            k=channel_intersections[[i]][j]
            print(c(i,k))
            hec_linestmp_old=hec_linestmp

            #@ Iteratively update hec_linestmp by inserting lateral weir
            #@ Use only unique parts of storage areas -- HEC RAS cannot have
            #@ overlapping lateral structures. Also, don't let lateral
            #@ structures span over bridges

            #@ TRICKS TO MAKE THINGS WELL BEHAVED.
            #@ Make sure that the chan_boundary_points we search through actually belong to this channel
            #@ , as sometimes we can get points intersecting multiple channels
            #@ 
            kk=which(!is.na(over(chan_boundary_points, chan2_list[[k]]))) # Indices of bank points within channel poly
            reachlist=sort( table( as.character(chan_boundary_points@data[kk,1])),decreasing=TRUE) # Names of reaches
            reachname=names(reachlist)[1] # Most commonly occurring name -- this will be the reach name
            chan_tmp_points=which(as.character(chan_boundary_points@data[kk,1])==reachname)
            chan_tmp_points=chan_boundary_points[kk[chan_tmp_points],]

            #@ We cannot have a lateral structure between the first and second points in a reach.
            #@ Also, we tend to have problems if we put it near the most downstream point.
            #@ So let's remove the first points from chan_tmp_points.
            remove=c(1, min(which(chan_tmp_points$bank=='R')))
            chan_tmp_points=chan_tmp_points[-remove,]
            remove=c( max(which(chan_tmp_points$bank=='L')), length(chan_tmp_points[,1]))
            chan_tmp_points=chan_tmp_points[-remove,]

            #@ We must ensure that the elevations along the lateral structure
            #@ are not below the min elevation in the storage area
            min_structure_elev=min(new_store_stage_vol_list[[i]][,1])
            print(c('#3322311# ', min_structure_elev))


            #@ 
            hec_linestmp=make_lateral_weir_text(unique_new_store_list[[i]], new_storage_names[[i]],
                                                     chan2_list[[k]], chan_tmp_points, 
                                                     lidar_DEM, vertical_datum_offset, 
                                                     hec_linestmp,
                                                     bridge_lines, 
                                                     limit_weir_elevation_by_channel_bank_elevation,
                                                     min_structure_elev, lower_limit_on_lateral_weir_elevations)

            # Trick to join up storage areas which are 'just' missing the channel point.
            l1=length(hec_linestmp_old)
            l2=length(hec_linestmp)
            l3=min(l1,l2)
            if(all(hec_linestmp_old[l3]==hec_linestmp[l3])& (l2==l1)){
                #browser()
                unique_new_store_list[[i]]=gBuffer(unique_new_store_list[[i]],width=20.) # Buffer by 20.
                hec_linestmp=make_lateral_weir_text(unique_new_store_list[[i]], new_storage_names[[i]],
                                                         chan2_list[[k]], chan_tmp_points, 
                                                         lidar_DEM, vertical_datum_offset, 
                                                         hec_linestmp,
                                                         bridge_lines, 
                                                         limit_weir_elevation_by_channel_bank_elevation,
                                                         min_structure_elev)


            }

        }

    }

    hec_lines2=hec_linestmp
    #@ Write to output
    cat(hec_lines2,file=output_file,sep="\n") 

    sink()
}
