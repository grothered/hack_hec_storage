# Utility functions for hec-help

pad_string<-function(mystring_vec, charlen, pad=" ", justify='left'){
    #@ FOR FORMATTING CHARACTER STRINGS TO HAVE A FIXED LENGTH
    #@ Take a string, check its length, and if it is too short, then pad it to
    #@ the correct length. If it is too long, then truncate it. 

    out=mystring_vec # Predefine
  
    # Loop over each vector element seperately 
    for(i in 1:length(mystring_vec)){ 
        mystring=mystring_vec[i]

        l=nchar(mystring) # Input string length

        if(l<=charlen){
            # String needs to be padded
            padtext=paste(rep(pad,charlen-l), collapse="")

            if(justify=='left'){
                out[i]=paste(mystring,padtext,sep="")
            }else if(justify=='right'){
                out[i]=paste(padtext,mystring,sep="")
            }

        }else if(l>charlen){
            # String needs to be truncated
            mystring_split=strsplit(mystring,"")[[1]] # Split to individual characters

            if(justify=='left'){
                out[i]=paste(mystring_split[1:charlen],sep="", collapse="")
            }else if(justify=='right'){
                out[i]=paste(mystring_split[(l-charlen+1):l], sep="", collapse="")
            }
        }
    }

    out
}

###########################################################

format_in_rows<-function(char_vec,rowlen){
    #@ Function to re-format a character vector to have rows containing 'rowlen' entries
    #@ e.g for rowlen=2
    #@ c(' 01', ' 23', ' 34.5', 'my char', 'f') --> c(' 01 23', ' 34.5my char', 'f') 

    l=length(char_vec)
    #@ l1= Number of lines we need to add to make the number of next_lines elements a multiple of rowlen
    l1= (rowlen-l%%rowlen)%%rowlen 

    char_vec=c(char_vec, rep(" ", l1)) # Padded

    #@ Convert to matrix with rowlen columns
    char_vec2=matrix(char_vec,ncol=rowlen,byrow=T)
    
    #@ Convert to rows with rowlen entries
    #e.g. char_vec3=paste(char_vec2[,1], char_vec2[,2], char_vec2[,3], char_vec2[,4], char_vec2[,5], sep="")
    char_vec3=c()
    for(i in 1:rowlen){ 
        char_vec3=paste(char_vec3, char_vec2[,i],sep="") 
    }

    char_vec3
}


split_nchars_numeric<-function(string,no_chars){
    #@ FOR READING STRUCTURED TEXT
    #@ Take a string, split it into individual characters, and recombine them
    #@ into separate pieces of no_chars length. Then convert these to numeric,
    #@ and output as a vector

    tmp=unlist(strsplit(string, split=""))
    l=length(tmp)
    if(l%%no_chars!=0){
        print(c('tmp = ', tmp))
        print(c('no_chars = ', no_chars))
        print(c('Length tmp = ', l))
        stop('ERROR in split_nchars_numeric: Number of characters in string is not a multiple of no_chars')

    }

    l2=l/no_chars
    output=rep(NA,l2)

    for(i in 1:l2){
        output[i] = as.numeric(paste(tmp[(1:no_chars) + no_chars*(i-1)], sep="", collapse=""))
    }
    output
}

########################################################################

create_channel_polygon<-function(hec_chan_file, spatial_proj){
    #@ Function to read a hec-ras file, extract the channel network, and make a
    #@ SpatialPolygonsDataFrame which covers it

    #library(sp)
    #library(rgdal)

    #@ Read file
    fin=file(hec_chan_file, open='r')
    hec_lines=readLines(fin)
    close(fin)

    #@ Identify channels
    chan_ind=grep('River Reach=', hec_lines)
    chan_polygons=list()

    #@ For each channel, extract the cross-sectional start/end points
    #@ lb = start (lower bound), ub = end (upper bound)
    for(i in 1:length(chan_ind)){
        offset=chan_ind[i]-1
        lb=chan_ind[i]
        if(i<length(chan_ind)){
            ub=chan_ind[i+1]
        }else{
            ub=length(hec_lines)
        }
        #@ Only grep the relevant cross-sections
        xsect_start=grep('XS GIS Cut Line=', hec_lines[lb:ub]) + offset
        #xsect_end=grep('Node Last Edited Time=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        xsect_end=grep('#Sta/Elev=', hec_lines[lb:ub]) + offset - 1

        # Treat the case of channels without xsections
        if(length(xsect_start)==0){
            err_mess = paste('ERROR: The following reach appears to not have any cross-sections', 
                             '\n', hec_lines[chan_ind[i]], '\n 
                               This will cause hec_help.R to fail :( \n',
                               ' Try adding 2 artificial cross-sections to the reach, or deleting it')
            stop(err_mess)
        }
        
        # Logical check that start and ends of x-sections line up.
        for(j in 1:length(xsect_start)){
            if(xsect_start[j]>xsect_end[j]){
                stop('ERROR: xsect_start is greater than xsect_end.  Check that
                      all cutlines in your hec-ras geometry file are georeferenced')
            }
        }
   
        #@ NOW EXTRACT THE START AND END POINTS OF EACH CROSS-SECTION 

        #@ Predefine variables to store start and end points of
        #@ cross-sections
        coords_start=c()
        coords_end=c()
        for(j in 1:length(xsect_start)){  # Loop over all xsections
            cutline_text=hec_lines[(xsect_start[j]+1):(xsect_end[j]-1)]
            coords=c()
            for(k in 1:length(cutline_text)){ # Loop over all lines 
                coords1=split_nchars_numeric(cutline_text[k],16)
                coords=rbind(coords,matrix(coords1,ncol=2,byrow=TRUE))
            }
            #print(coords)
            coords_start=rbind(coords_start,coords[1,])
            coords_end=rbind(coords_end,coords[length(coords[,1]),])
        }

        #@ Combine coordinates and make a closed polygon
        coords_all=rbind(coords_start, coords_end[length(coords_end[,1]):1,], coords_start[1,])

        #@ Append to list of polygons, in a way that makes it easy to produce a
        #@ spatialpolygonsdataframe later

        chan_polygons[[i]] = (Polygons(list(Polygon(coords_all)), ID=as.character(i)))
        if(chan_polygons[[i]]@area ==0){
            print( 
                 paste(c( '##############################################################################\n',
                          'ERROR: zero area for channel polygon ', i,' , ', hec_lines[chan_ind[i]], '\n',
                          'This probably means that the reach has less than 2 cross-sections \n',
                          'You need at least 2 cross-sections on every reach \n',
                          '##############################################################################\n',
                          ' '))
                 )
            chan_polygons[[i]]=NULL
        }
        
    }
        # Catch errors in reach definition
        if(any( unlist ( lapply(chan_polygons, is.null) ) ) ){
            stop('Will not proceed until you fix the above reaches')
        }
        #@ Now coerce output to SpatialPolygonsDataFrame
        chan=SpatialPolygons(chan_polygons, proj4string=spatial_proj)

        chan2=SpatialPolygonsDataFrame(chan, data=data.frame(id=seq(1,length(chan_polygons))), match.ID=FALSE)


        chan2
}

#################################################################################################
get_existing_storage_areas<-function(hec_file,spatial_proj){
    #@ Function to extract storage areas from an existing hec-ras file
    #@ Useful to plot them up, to help with the creation of new storage areas.
    fin=file(hec_file, open='r')
    hec_lines=readLines(fin)
    close(fin)
   
    store_coords_start=grep('Storage Area Surface Line', hec_lines) +1
    store_coords_end=grep('Storage Area Type', hec_lines) -1

    #@ Treat case with no storage areas
    if(length(store_coords_start)==0) return(NULL) 

    #@ Loop over all storage areas, and add them to a SpatialPolygonsDataFrame
    poly_list=list()
    storage_names=c()
    for(i in 1:length(store_coords_start)){
        storage_name=hec_lines[store_coords_start[i]-2]
        storage_name=strsplit(storage_name, '=')[[1]][2]
        storage_name=strsplit(storage_name,',')[[1]][1]
        
        coord_text=hec_lines[store_coords_start[i]:store_coords_end[i]]
        coords_out=c()
        for(j in 1:length(coord_text)){
            #print(coord_text[j])
            coords_out=c(coords_out,split_nchars_numeric(coord_text[j],16))
        }
        coords_out=matrix(coords_out,ncol=3,byrow=T)
        coords_out=rbind(coords_out, coords_out[1,])

        poly_list[[i]]=Polygons(list(Polygon(coords_out[,1:2])), ID=storage_name)
        storage_names=c(storage_names, storage_name) 
    }
    storage_sp_poly=SpatialPolygons(poly_list, proj4string=CRS(spatial_proj))

    storage_sp_polydf=SpatialPolygonsDataFrame(storage_sp_poly, data=data.frame(name=storage_names), match.ID=FALSE)
    storage_sp_polydf
}
#################################################################################################

compute_stage_vol_relation<-function(my_poly,lidar_DEM, vertical_datum_offset, upper_bound_stage=100){
    #@
    #@ Function to compute the stage-volume relation for a region defined by the polygon my_poly
    #@ , on a given lidar_Dem raster  
    #@ Output: A Stage-Volume relation, with the stage values offset by + vertical_datum_offset
    #@

    cell_area=xres(lidar_DEM)*yres(lidar_DEM)/1000. # Area of a pixel -- units should match the units of the output data.
    #@ In hec-ras, the storage area units are in 1000s of m^3 -- so we need 10
    #@ lots of 10*10 squares filled in a metre of water to account for one unit.
    #@ Hence the strange use of units above

    #elev_pts=lidar_DEM[my_poly] # Elevation points inside the polygon
    elev_pts=extract(lidar_DEM, my_poly, small=TRUE)
   
    #@ Number of points on the stage - volume curve. 
    #@ Ensure 5<= number of points <= 60. Ideally, only give a new point every 0.33 m 
    num_stagevol_points=min(30, max(5, (max(elev_pts[[1]])-min(elev_pts[[1]]) )*3 ) ) 

    elev_hist=hist(elev_pts[[1]],n=num_stagevol_points)

    # Compute a sequence of stages at which we will evaluate the stored volume
    lower_bound=elev_hist$mids[1]-0.1 #
    upper_bound=max(upper_bound_stage,max(elev_hist$mids)+0.01) # MAX value at which we want the volume. This is a trick so that hec-ras never exceeds the range
    desired_stage=c(lower_bound,elev_hist$mids, upper_bound) # This is the range of stages for which we want output
    area_fun=approxfun(desired_stage, c(0,cumsum(elev_hist$counts)*cell_area, sum(elev_hist$counts)*cell_area), rule=2 )
    # Function to compute volume (note: volume = integral (Area))
    vol_fun<-function(stage){
        integral=integrate(area_fun, lower_bound, stage,subdivisions=1000)
        integral$value
    }

    output=cbind(desired_stage+vertical_datum_offset, sapply(desired_stage,vol_fun))

    output
}

###################################################################################################

make_storage_area_text<-function(storage_area, elev_vol,name){
    #@ Function to make a block of text representing a storage area,
    #@ as appears in the hec-ras format .g01
    #@
    #@ e.g.
    #@
    #@ Storage Area=east_bank1      ,510988.0387824,1620987.9181366
    #@ Storage Area Surface Line= 15 
    #@  510643.523093591621316.98005342                
    #@ 510548.0930744711621305.29474495                
    #@  510709.739841551621272.18637097                
    #@ 510723.3727014251621141.70042646                
    #@ 510738.9531127091620997.58162208                
    #@ 510735.0580098881620913.83691142                
    #@ 510700.0020844971620832.03975218                
    #@ 510672.7363647491620759.98034998                
    #@ 510657.1559534641620682.07829356                
    #@ 511175.2046286831620545.74969482                
    #@ 511336.8513957621620814.51178948                
    #@ 511453.7044803981621168.96614621                
    #@ 511299.8479189611621264.39616533                
    #@ 510898.6523283781621262.44861392                
    #@ 510666.8937105171621313.08495059                
    #@ Storage Area Type= 1 
    #@ Storage Area Area=
    #@ Storage Area Min Elev=
    #@ Storage Area Vol Elev= 84 
    #@    17.35       0   17.45    .005   17.55    .025   17.65    .075   17.75    .185
    #@    17.85.3600001   17.95.6250002   18.051.019999   18.15   1.595   18.252.375064
    #@    18.353.335002   18.454.444993   18.555.669004   18.656.999982   18.758.465041
    #@    18.8510.15023   18.95   12.17   19.0514.55502   19.1517.26011   19.2520.21989
    #@    19.3523.39999   19.4526.79979   19.5530.42019   19.6534.26412   19.75  38.299
    #@    19.85   42.64   19.9547.40543   20.0552.62051   20.1558.45985   20.25 65.0401
    #@    20.3572.31497   20.4580.16493   20.5588.45006   20.6597.16591   20.75106.3653
    #@    20.85116.1453   20.95126.5704   21.05137.5099   21.15 148.912   21.25160.7833
    #@    21.35173.0752   21.45 185.787   21.55198.9568   21.65212.6048   21.75226.5446
    #@    21.85240.7119   21.95255.1156   22.05269.7048   22.15 284.454   22.25299.3548
    #@    22.35314.3996   22.45329.6109   22.55344.9856   22.65360.5155   22.75376.2579
    #@    22.85392.1823   22.95408.4679   23.05424.5935   23.15 441.215   23.25458.2499
    #@    23.35475.6863   23.45493.6751   23.55512.3254   23.65532.1439   23.75553.4506
    #@    23.85577.0807   23.95604.1621   24.05634.8083   24.15668.5247   24.25 704.508
    #@    24.35741.9367   24.45780.8355   24.55 820.769   24.65861.5136   24.75904.3318
    #@    24.85944.7191   24.95986.9531   25.051029.436   25.151072.018   25.251114.711
    #@    25.351157.568   25.451200.513   25.551243.467   110.537729.49
    #@
    #@
    #@ 
    #@    These string manipulations can be pretty hacky
    #@

    # Try to filter storage area coordinates -- we sometimes have way too many points
    storage_area=gSimplify(storage_area,tol=10., topologyPreserve=TRUE)

    #@ Compute coordinates inside poly
    storage_area_central_coords=coordinates(storage_area)

    #@ Compute bounding coordinates -- got to love the notation
    storage_area_bounding_coords=coordinates(gBoundary(storage_area))[1][[1]][[1]]

    #@ Remove final bounding coordinate (which is a repeat of the first point, and does not feature in hec)
    l=length(storage_area_bounding_coords[,1])
    storage_area_bounding_coords=storage_area_bounding_coords[-l,]

    ## TRY TO FILTER storage_area_bounding_coords as stored in hec-ras
    ## Because it can be too large often, perhaps due to the details of the spatial processing
    ## Hopefully not needed, as we simplify above
    if(length(storage_area_bounding_coords[,1])>50.){
        print('More than 50 points still in storage area')
        usdist=c(0, cumsum(sqrt(diff(storage_area_bounding_coords[,1])**2 + diff(storage_area_bounding_coords[,2])**2)))
        x_fun=approxfun(usdist,storage_area_bounding_coords[,1])
        y_fun=approxfun(usdist,storage_area_bounding_coords[,2])

        lmax=max(usdist)
        # How many coordinates can we accept? Max of 500. Min 10m spacing. 
        num_coords=min(499, floor(lmax/10.), length(storage_area_bounding_coords[,1]))
        dist_seq=seq(0,lmax,len=num_coords)
        x_new=x_fun(dist_seq)
        y_new=y_fun(dist_seq)

        storage_area_bounding_coords=cbind(x_new,y_new)
    }

    #browser()
    #@ Start producing output text

    output_text=c() # Initialise output_text

    #@ Create character strings of correct format for name
    name_nonempty=as.character(name)
    name_fil=rep(" ", 16-nchar(name_nonempty))
    name_fil=paste(name_fil,collapse='')

    #@ Create coordinate character strings with correct format
    coord1=format(as.character(signif(storage_area_central_coords[1,1],13)),width=15,justify='none',trim=T)
    coord2=format(as.character(signif(storage_area_central_coords[1,2],13)),width=15,justify='none',trim=T)

    #@ Build first line of output text
    first_line=paste("Storage Area=", name_nonempty, name_fil, ",", coord1,",", coord2, sep="")
    #print(first_line)
    output_text=c(output_text, first_line)

    #@ Build the second line of output text
    second_line=paste('Storage Area Surface Line=', l-1)
    output_text=c(output_text,second_line)

    #@ Build boundary output coords with correct format
    #coord1=format(as.character(signif(storage_area_bounding_coords[,1],16)),width=17,justify='none',trim=T)
    #coord2=format(as.character(signif(storage_area_bounding_coords[,2],16)),width=17,justify='none',trim=T)
 
    coord1=pad_string(as.character(signif(storage_area_bounding_coords[,1],15)), 16,pad="0",justify='left') 
    coord2=pad_string(as.character(signif(storage_area_bounding_coords[,2],15)), 16,pad="0",justify='left') 
 
    #@ Forcibly pad these - it is critical to get the width right
    #FIXME: Must be a more elegant way
    #n1=nchar(coord1)
    #for(i in 1:length(n1)){
    #    if(n1[i]<16){
    #        coord1[i]=paste(paste(rep(" ", 16-n1[i]),collapse=""), coord1[i],sep="")
    #    }
    #}
    
    next_lines=paste(coord1,coord2,sep="")
    l=length(coord1)
    next_lines=paste(next_lines, rep(" ", l, collapse=""),sep="") 
    output_text=c(output_text,next_lines)

    #@ Add storage area type info
    output_text=c(output_text, c("Storage Area Type= 1", "Storage Area Area=","Storage Area Min Elev="))

    #@ Add volume - elev information
    l=length(elev_vol[,1])
    output_text=c(output_text,paste('Storage Area Vol Elev=', l) )
    #elev_coord=format(as.character(round(elev_vol[,1],2)), width=8,justify='right',trim=T)
    elev_coord=pad_string(as.character(round(elev_vol[,1],2)), charlen=8,pad=" ", justify='right')
    #vol_coord=elev_vol[,2]
    #vol_coord=format(as.character(signif(elev_vol[,2],5)),width=8,justify='right',trim=T)
    vol_coord=pad_string(as.character(signif(elev_vol[,2],6)), charlen=8,pad=" ", justify='right')

    next_lines=paste(elev_coord,vol_coord,sep="")
    # Get the format right
    next_lines = format_in_rows(next_lines,5)

    output_text=c(output_text,next_lines)

    output_text
    #print(output_text)
}

####################################################################################################################

make_storage_connection_text<-function(store1, store2, name1, name2, lidar_DEM, vertical_datum_offset=10.5){
    #@
    #@ Function to develop the text for a storage area connection between 2 overlapping polygons, store1 & store2
    #@ It finds the hypsometry in their overlapping zone, and uses this to compute the geometry of their weir
    #@ The text is for insertion into a hec-ras .g01 file
    #@

    #@ The format looks like this
    #@
    #@ Connection=south_east_mont_,514493.3662094,1627706.738947
    #@ Connection Desc=south_east_mont_1_to_2
    #@ Connection Line= 3 
    #@ 514532.9213419731627758.28048332514532.9213419731627758.28048332
    #@ 514414.2559442781627603.65587421
    #@ Connection Last Edited Time=Jun/03/2012 19:27:42
    #@ Connection Up SA=south_east_mont1
    #@ Connection Dn SA=south_east_mont2
    #@ Conn Routing Type= 1 
    #@ Conn Use RC Family=True
    #@ Conn Weir WD=20
    #@ Conn Weir Coef=1.66
    #@ Conn Weir Is Ogee= 0 
    #@ Conn Simple Spill Pos Coef=0.05
    #@ Conn Simple Spill Neg Coef=0.05
    #@ Conn Weir SE= 2 
    #@        0      24     300      24
    #@ Conn HTab HWMax=40
    #@
    #@
    print(' ')
    print('#########################################')
    print(paste('Connecting ', name1, 'with', name2))
    print(paste('These are located near'))
    print(coordinates(store1))
    print(coordinates(store2))
    print('#########################################')
    print(' ')    
    intersection=gIntersection(store1, store2) # Polygon containing the intersection of the 2 storage areas
   
    #@ Check that they do not overlap too much -- could be a problem
    A1=gArea(store1)
    A2=gArea(store2)
    A3=gArea(intersection)
    
    overlap_max=A3/min(A1,A2)
    overlap_error_threshold=0.2
    if(overlap_max>overlap_error_threshold){
        if(A1>A2){
            plot(store1,axes=T,asp=1)
            plot(store2,border='red',add=T)
            title('Storage areas which overlap too much')
        }else{
            plot(store2,axes=T,asp=1)
            plot(store1,border='red',add=T)
            title('Storage areas which overlap too much')
        }
        print('')
        print('#################################################')
        stop(paste('ERROR: Storage areas', name1, 'and', name2, 'overlap by > ', 
                    overlap_error_threshold*100, '%. This sounds like an error', 
                    '- check your input data'))
    }
    
    output_text=c() # Predefine output

    #@
    #@ Compute coordinates for connection
    #@
    connection_pts=coordinates(intersection)
    #@ Format connection_pts
    con_x=pad_string(as.character(connection_pts[1]), 16,pad="0")
    con_y=pad_string(as.character(connection_pts[2]), 16,pad="0")
    
    #@
    #@ Make name for connection -- include a timestamp to make unique
    #@
    name_timestamp=as.character(floor(as.numeric(Sys.time())*100)%%100000000)
    #connection_name=paste('SC______',name_timestamp,sep="") 
    connection_name=paste(name1,name2,sep="")
    long_connection_name=paste('Connect', name1, '&', name2)

    first_line=paste('Connection=', connection_name,",", con_x,",", con_y,sep="")
    second_line=paste('Connection Desc=', long_connection_name, sep="")
    third_line= paste('Connection Line= 2')

    #@
    #@ now make up coordinates for connection line
    #@

    #@ Line direction = direction of line joining the 2 centroids
    store1_cent=coordinates(store1)
    store2_cent=coordinates(store2)
    vec=store1_cent-store2_cent
    vec=vec/10 # This vector will give us the length of the drawn storage area connection line
    p1=connection_pts+vec
    p2=connection_pts-vec
    #@ Format line points
    p1x=pad_string(as.character(p1[1]), 16, pad="0")
    p1y=pad_string(as.character(p1[2]), 16, pad="0")
    p2x=pad_string(as.character(p2[1]), 16, pad="0")
    p2y=pad_string(as.character(p2[2]), 16, pad="0")
    fourth_line=paste(p1x, p1y, p2x, p2y, sep="")
    
    #@
    #@ Now make the 'connection_edited' timestamp
    #@
    edit_time=format(Sys.time(), "%b/%d/%Y %H:%M:%S") # Timestamp like hec-ras
    fifth_line=paste("Connection Last Edited Time=", edit_time,sep="")

    output_text=c(first_line, second_line, third_line, fourth_line, fifth_line)

    #@
    #@ Make connection info
    #@
    next_lines=c(paste('Connection Up SA=', name1, sep=""), 
                 paste('Connection Dn SA=', name2, sep=""))

    output_text=c(output_text,next_lines)
  
    #@ 
    #@ Add other weir info
    #@ 
    next_lines=c("Conn Routing Type= 1",  
                 "Conn Use RC Family=True",
                 "Conn Weir WD=20",
                 "Conn Weir Coef=1.66",
                 "Conn Weir Is Ogee= 0", 
                 "Conn Simple Spill Pos Coef=0.05",
                 "Conn Simple Spill Neg Coef=0.05")
    
    output_text=c(output_text,next_lines) 


    #@
    #@ Compute weir form 
    #@
    
    #elev_relation=lidar_DEM[intersection][[1]]
    elev_relation=extract(lidar_DEM, intersection, small=TRUE)[[1]]
    #if(length(elev_relation)==0){
    #    print(' ')
    #    print('WARNING: Not creating a storage area connection for')
    #    print(paste(name1, ' to ', name2))
    #    print('Because the overlap is too small to enclose any points in the DEM')
    #    return(NA)
    #}
        
    elev_relation=elev_relation+ vertical_datum_offset
    elev_relation2=sort(elev_relation)
    #@ How to compute the weir length??
    #@ Idea:
    #@ gLength(intersection) = poly boundary length
    #@ If we assume poly is long and thin, then length~ = boundary length /2 
    #@ We could make it of slightly shorter length to be conservative (e.g. boundary_length *1/3 or *5/12)
    weir_length = gLength(intersection)*1/3
   
    minimum_weir_length=5 
    if(weir_length<minimum_weir_length){
        print(' ')
        print('WARNING: Not creating a storage area connection for')
        print(paste(name1, ' to ', name2))
        print('Because the computed weir length is small. It is')
        print(weir_length)
        return(NA)
    }
        
    l=length(elev_relation)
    weir_x_vals=seq(0,weir_length, len=l) # X values at which we will get weir elevation points
    l2=min(l,30) # Number of points on the weir in hec-ras

    if(l>1){
        weir_relation=approx(weir_x_vals, elev_relation2, n=l2) # Weir x - elev relation
    }else{
        # Treat the case of only one point
        l2=2
        weir_relation=approx(c(0, weir_length), c(elev_relation2, elev_relation2),n=l2)
    }

    weir_x=pad_string(as.character(signif(weir_relation$x,7)), 8, pad=" ", justify='right')
    weir_y=pad_string(as.character(signif(weir_relation$y,7)), 8, pad=" ", justify='right')

    weir_text=paste(weir_x,weir_y,sep="")
    weir_text=format_in_rows(weir_text,5)
    output_text=c(output_text, paste("Conn Weir SE=", l2))
    output_text=c(output_text, weir_text)


    #@ Append the Htab HWMax parameter
    output_text=c(output_text, paste("Conn HTab HWMax=", round(max(elev_relation2)+15,2),sep=""))

    output_text
}

###############################################################################################################

make_channel_boundary_points<-function(hec_chan_file,spatial_proj){
    #@ Produce a shapefile containing the end points of the cross-sections 

    #@ Read file
    fin=file(hec_chan_file, open='r')
    hec_lines=readLines(fin)
    close(fin)

    #@ Identify channels
    chan_ind=grep('River Reach=', hec_lines)

    chan_polygons=list()

    #@ For each channel, extract the cross-sectional start/end points
    output_coords=c()
    output_data=c()
    for(i in 1:length(chan_ind)){

        offset=chan_ind[i]-1
        lb=chan_ind[i]
        if(i<length(chan_ind)){
            ub=chan_ind[i+1]
        }else{
            ub=length(hec_lines)
        }

        #@ Only grep the relevant cross-sections
        xsect_start=grep('XS GIS Cut Line=', hec_lines[lb:ub]) + offset
        #xsect_end=grep('Node Last Edited Time=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        xsect_end=grep('#Sta/Elev=', hec_lines[lb:ub]) + offset - 1
      
        if(length(xsect_start)==0){
            err_mess = paste('ERROR: The following reach appears to not have any cross-sections', 
                             '\n', hec_lines[chan_ind[i]], '\n 
                               This will cause hec_help.R to fail :( \n',
                               ' Try adding 2 artificial cross-sections to the reach, or deleting it')
            stop(err_mess)
        }
        # Useful to have the 'station' lines as well. Initially I tried to use
        # this as the definition of xsect_start, however, I decided that was
        # not the best. However, we need this line to get the station number
        station_start=grep('Type RM Length L Ch R =', hec_lines[lb:ub]) + offset
   
        #@ NOW EXTRACT THE START AND END POINTS OF EACH CROSS-SECTION 

        #@ Predefine variables to store start and end points of
        #@ cross-sections
        coords_start=c()
        coords_end=c()
        reach_name=c()
        station_name=c()
        left_bank_downstream_dist=c()
        right_bank_downstream_dist=c()
        left_bank_elev=c()
        right_bank_elev=c()
        line_num=c()
        #@ Loop over all xsections and extract required data
        for(j in 1:length(xsect_start)){  
            cutline_text=hec_lines[(xsect_start[j]+1):(xsect_end[j]-1)]
            #@ Get cutline coordinates
            coords=c()
            for(k in 1:length(cutline_text)){ # Loop over all lines 
                coords1=split_nchars_numeric(cutline_text[k],16)
                coords=rbind(coords,matrix(coords1,ncol=2,byrow=TRUE))
            }
            coords_start=rbind(coords_start,coords[1,])
            coords_end=rbind(coords_end,coords[length(coords[,1]),])
            # Get reach name
            reach_name=c(reach_name, strsplit(hec_lines[offset+1], "=")[[1]][2] )
            # Get station name
            station_index=max(station_start[station_start<=xsect_start[j]])
            station_name=c(station_name, strsplit(hec_lines[station_index], ",")[[1]][2] )
            
            # Get left and right bank downstream distances
            left_bank_downstream_dist=c(left_bank_downstream_dist, strsplit(hec_lines[station_index], ",")[[1]][3] )
            right_bank_downstream_dist=c(right_bank_downstream_dist, strsplit(hec_lines[station_index], ",")[[1]][5] )
            
            # Get elevation associated with left and right banks
            sei=xsect_end[j]+1 # Index of the '#Sta/Elev' line, just before the station-elevation data
            # FIXME: Here we are assuming that there are < 300 station elevation lines -- probably true, but be careful!
            mi = grep('#Mann=',hec_lines[sei:(sei+300)])[1]+sei-1 # Index of the #Mann= line after the station-elevation data 
            xsect_sta_elev_inds=(sei+1):(mi-1) # Indices of the station - elevation data
            xsect_station_elevation=split_nchars_numeric(hec_lines[xsect_sta_elev_inds], 8.) # Each number takes up 8 characters
            left_bank_elev = c(left_bank_elev, xsect_station_elevation[2])
            right_bank_elev = c(right_bank_elev, xsect_station_elevation[length(xsect_station_elevation)])
            line_num = c(line_num,xsect_start[j])
        }
       
        output_coords=rbind(output_coords, coords_start)
        output_data=rbind(output_data, 
                          cbind(reach_name, station_name, 
                                rep('L', length(reach_name)), 
                                left_bank_downstream_dist,
                                left_bank_elev,
                                line_num))
        output_coords=rbind(output_coords, coords_end)
        output_data=rbind(output_data, cbind(reach_name, station_name, 
                                             rep('R', length(reach_name)), 
                                             right_bank_downstream_dist, 
                                             right_bank_elev,
                                             line_num))
        
        #print(cbind(station_name, reach_name))
    
    }
        # Coerce to spatial points
        output_pts=SpatialPointsDataFrame(coords=output_coords[,1:2], 
                      data=data.frame(reach_name=output_data[,1], station_name=output_data[,2],
                                      bank=output_data[,3], downstream_distance=as.numeric(output_data[,4]),
                                      bank_elev=as.numeric(output_data[,5]), line_num=as.numeric(output_data[,6])),
                      match.ID=FALSE,
                      proj4string=CRS(spatial_proj))
        output_pts
}

##############################################################################################################
make_channel_cutlines<-function(hec_chan_file, spatial_proj){
    #@ Make a spatialLines object containing the cross-sectional cutlines

    #@ Read input file
    fin=file(hec_chan_file, open='r')
    hec_lines=readLines(fin)
    close(fin)

    #@ Identify cutlines
    cutline_start=grep('XS GIS Cut Line', hec_lines)+1
    cutline_end=grep('#Sta', hec_lines)-2

    lines_list=list()
    for(i in 1:length(cutline_start)){
        cutline=split_nchars_numeric(hec_lines[cutline_start[i]:cutline_end[i]], 16)
        cutline_2=matrix(cutline,byrow=T,ncol=2)

        lines_list[[i]] = Lines(list(Line(cutline_2)), ID=as.character(i))
    }

    #@ Identify associated reaches
    reaches=grep('River Reach', hec_lines)
    xsect_reaches=cutline_start*NA
    for(i in 1:length(xsect_reaches)){
        xsect_reaches[i] = reaches[ which.max( cumsum(reaches< cutline_start[i]) ) ]
    }

    xsect_labels=grep('Type RM Length L Ch R = 1', hec_lines) 

    xsect_cutlines=SpatialLines(lines_list,proj4string=CRS(spatial_proj))
    output=SpatialLinesDataFrame(xsect_cutlines, data=data.frame(id_2=1:length(xsect_cutlines), 
                                 sec_info=as.character(hec_lines[xsect_labels]), 
                                 reach=as.character(hec_lines[xsect_reaches])), match.ID=FALSE)
    return(output)
}

###########################################################################################################################################

compute_centrelines<-function(hec_chan_file,spatial_proj){
    #@ Extract channel centrelines from hecras file

    #@ Read input file
    fin=file(hec_chan_file, open='r')
    hec_lines=readLines(fin)
    close(fin)

    #@ Identify reaches
    reaches=grep('River Reach', hec_lines)
    reach_coordinate_end=grep('Rch Text', hec_lines)

    if(length(reaches)!=length(reach_coordinate_end)){
        stop("ERROR in compute_new_downstream_distances: The start and end points of the reach coordinates are not matching up")
    }

    #@ Now extract the centreline coordiantes
    centrelines=list()
    for(i in 1:length(reaches)){
        centre_txt=hec_lines[(reaches[i]+2):(reach_coordinate_end[i]-1)]
        centre_coords=matrix(split_nchars_numeric(centre_txt, 16), ncol=2,byrow=T)    

        centrelines[[i]] = Lines(list(Line(centre_coords)), ID=hec_lines[reaches[i]])
    }
    
    reach_lines=SpatialLines(centrelines,proj4string=CRS(spatial_proj))
    reach_lines=SpatialLinesDataFrame(reach_lines, data=data.frame(id_2=1:length(reach_lines), reach=as.character(hec_lines[reaches])), match.ID=FALSE)
    return(reach_lines) 
}

#########################################################################################################################################

update_downstream_distances<-function(hec_lines, chan_cutlines, reach_lines){
    #@ Function to compute the 'downstream distances' of the xsections on each
    #@ reach, using the cutline and reach centreline information in hec_lines.
    
    #@ Assume that left bank = right bank = straight line 
    #@ distance between x-sections at intersection with the channel,
    #@ while channel distance = distance along the channel.

    # NOTE: We are setting the left bank and right bank distances to be the same.
    # This should be appropriate if the overbank flows are largely occurring near the channel
    # However, this will not always be correct, although it will often be true
    # if the floodplain roughness gets high away from the channel. 
    # In reality, the left bank / right bank distances should depend on the
    # flow state -- during serious inundation, they might differ quite a bit. 
    
    spatial_proj=proj4string(chan_cutlines)

    # Now, for every reach_line, find its intersection with the centreline.
    # This should be just one point: if more than that, we need to do some more
    # work
    #downstream_distances_straight=list()
    #downstream_distances_chan=list()

    for(i in 1:length(reach_lines)){
        # Get cross-sections located on this reach
        local_xsects_inds=which(chan_cutlines@data$reach==reach_lines@data$reach[i])
        local_xsects = SpatialLines(chan_cutlines@lines[local_xsects_inds], proj4string=CRS(spatial_proj))
        local_xsects_df = SpatialLinesDataFrame(local_xsects, data=chan_cutlines@data[local_xsects_inds,])#, proj4string=CRS(spatial_proj))
      
        # Get the reach line as a SpatialLines 
        local_centreline=SpatialLines(reach_lines@lines[i], proj4string=CRS(spatial_proj))

        #cutpoints=gIntersection(local_centreline,local_xsects,byid=T)

        # Loop over the xsections, and find their upstream_distance with the
        # channel, and intersection_point
        intersect_pt=c()
        upstream_distances_chan=rep(NA, length(local_xsects))
        for(j in 1:length(local_xsects)){
            cutpointz=gIntersection(local_xsects[j], local_centreline)

            # Treat the case of multiple intersections by choosing the one that
            # is most upstream
            if(length(cutpointz)>1){
                # Compute upstream distances of points along channel
                usdists=rep(NA,length(cutpointz))            
                for(k in 1:length(usdists)){
                    tmp = usdistfun(coordinates(cutpointz[k]), coordinates(local_centreline)[[1]][[1]])
                    usdists[k]=tmp[1]
                }
                # Select the point with the largest upstream distance
                keep_ind=which.max(usdists)
                cutpointz=cutpointz[keep_ind]
                upstream_distances_chan[j]= usdists[keep_ind]
                
            }else if(length(cutpointz)==0){
                print('no cutpoint')
                browser()
                stop('no cutpoint')
            }else{
                tmp = usdistfun(coordinates(cutpointz), coordinates(local_centreline)[[1]][[1]])
                upstream_distances_chan[j]=tmp[1] 
            }
            intersect_pt=rbind(intersect_pt, coordinates(cutpointz))
        }
        
        downstream_distances_chan=c(-diff(upstream_distances_chan), 0.)
        downstream_distances_straight=c((diff(intersect_pt[,1])**2 + diff(intersect_pt[,2])**2)**0.5, 0)
      
        # Now correct the hec_lines file 
        reach_index=grep(as.character(reach_lines@data$reach[i]), hec_lines)
        l = length(hec_lines)
        for(j in 1:length(local_xsects)){
            old_xsect_txt=as.character(local_xsects_df@data[j,]$sec_info) # Text string corresponding to the cross_section
            index_in_hecfile=grep(old_xsect_txt,hec_lines[reach_index:l], fixed=TRUE)[1] + reach_index-1 # Ensure that it occurs within the reach xsections

            text_split=strsplit(old_xsect_txt,",")[[1]]

            # Change the downstream distances
            if(downstream_distances_chan[j]!=0.){
                text_split[c(3,5)]=as.character(round(downstream_distances_straight[j],2))
                chan_dist=max(downstream_distances_straight[j], downstream_distances_chan[j])
                text_split[4] = as.character(round(chan_dist,2))
                

                text_combine=paste(text_split,collapse=",")
                hec_lines[index_in_hecfile]=text_combine
            }
                        
        }
    }

    return(hec_lines)
}

##############################################################

make_lateral_weir_text<-function(storage_poly, storage_name, chan_poly, chan_pts, lidar_DEM, 
                                 vertical_datum_offset=10.5, hec_lines, bridge_lines,
                                 limit_weir_elevation_by_channel_bank_elevation,
                                 min_structure_elev, lower_limit_on_lateral_weir_elevations){
    #@ Given a storage polygon and a channel polygon (which intersect), make me
    #@ a lateral weir and place it into hec_lines
    #@
    #@ Lateral weirs look like this:
    #@
    #@ Type RM Length L Ch R = 6 ,26250   ,,,
    #@ BEGIN DESCRIPTION:
    #@ south_east_mont2_to_channel
    #@ END DESCRIPTION:
    #@ Node Last Edited Time=May/04/2012 13:57:15
    #@ Lateral Weir Pos= 0 
    #@ Lateral Weir End=                ,                ,        ,south_east_mont2
    #@ Lateral Weir Distance=10
    #@ Lateral Weir TW Multiple XS=-1
    #@ Lateral Weir WD=20
    #@ Lateral Weir Coef=1.1
    #@ Lateral Weir WSCriteria=-1 
    #@ Lateral Weir Flap Gates= 0 
    #@ Lateral Weir Hagers EQN= 0 ,,,,,
    #@ Lateral Weir SS=0.05,0.05,
    #@ Lateral Weir Type= 0 
    #@ Lateral Weir Connection Pos and Dist= 0 ,
    #@ Lateral Weir SE= 2 
    #@        0    24.5    1050    24.5
    #@ Lateral Weir HW RS Station=26275.4*,-10
    #@ Lateral Weir TW RS Station=,0
    #@ LW Div RC= 0 ,False,
    #@
    #@

    #@ Compute intersection of storage poly with chan poly
    intersection=gIntersection(storage_poly,chan_poly)
    
    #@ Make sure that intersection contains only 1 connected polygon
    if(length(intersection@polygons[[1]]@Polygons)>1){
       #@ The intersection consists of > 1 disjoint polygon. Choose the one with largest area
       areas=lapply(intersection@polygons[[1]]@Polygons, getarea<-function(x) x@area) 
       areas=unlist(areas)
       tmp=which.max(areas) # Index of one with largest area
       tmp_Polygons=Polygons(list(intersection@polygons[[1]]@Polygons[[tmp]]), ID=1)
       intersection=SpatialPolygons(list(tmp_Polygons), proj4string=CRS(proj4string(storage_poly)))
    }
  
    #@ Get channel boundary points 
    vv=over(chan_pts, intersection)
    vv=which(!is.na(vv))
    if(length(vv)==0){
        print('##')
        print(paste('WARNING: There is an intersection of', storage_name, ' with the channel.'))
        print('However, it does not contain any bank points')
        print('Skipping this one')    
        print(' ')
        return(hec_lines)
    }
    weir_pts=chan_pts[intersection,]

    #@ If there are no points, return
    if(length(weir_pts)==0){
        return(hec_lines)
    }

    #@ Check to ensure we intersect the channel on only 1 bank
    banks=weir_pts$bank
    bank2=union(banks,banks)
    if(length(bank2)>1){
        print('####################')
        print(paste('WARNING: Storage polygon,', storage_name, "intersects on both the left and right banks"))
        print(paste('It is located near'))
        print(coordinates(weir_pts[1,]))
        num_left=sum(weir_pts$bank=='L')
        num_right=sum(weir_pts$bank=='R')
        print(paste('There are', num_left, 'points on the left bank, and', num_right, 'points on the right bank'))
        if(num_left>num_right){
            print('Using only points on the left bank')
            weir_pts=weir_pts[weir_pts$bank=='L',]
            bank2='L'
        }else{
            print('Using only points on the right bank')
            weir_pts=weir_pts[weir_pts$bank=='R',]
            bank2='R'
        }
        print(' ')
    }
    
    ##@ Check to ensure that the weir_pts do not include a bridge.
    ##@ Hecras cannot have lateral structure spanning a bridge
    weir_lines=weir_pts@data$line_num # Line numbers associated with the weir_lines xsections
    upstream_bridges=which(bridge_lines<min(weir_lines))
    downstream_bridges=which(bridge_lines> max(weir_lines))
    if(length(upstream_bridges)+length(downstream_bridges) < length(bridge_lines)){
        # Check that the weir points are sorted (they should be)
        if(min(diff(weir_lines))<0){
            stop('ERROR: Found unsorted weir_pts which apparently intersect a bridge. Need to code this case')
        }
        # Categorise weir_pts into groups without internal bridges
        weir_cut = cut(weir_lines, bridge_lines) 
        weir_sublen= sort(table(weir_cut), decreasing=TRUE)
        # Extract the most commonly occuring points
        longest_weir=which(weir_cut==as.factor(names(weir_sublen))[1])
        # Update the weir pts
        weir_pts=weir_pts[longest_weir,]
        lwp=length(weir_pts)
        if(lwp>2) weir_pts=weir_pts[2:(lwp-1),] # 
    }
    
    # Don't allow > 100 intersecting channel points
    if(length(weir_pts[,1])>99){
        weir_pts=weir_pts[1:99,]
    }
    
    
    #@ Sort the bank stations
    st1=as.character(weir_pts$station_name)
    st2=gsub('\\*', '', st1) # Remove * symbol
    st3=as.numeric(st2) # Now they are numbers

    #@ Sort the weir points along the channel
    station_order=sort(st3,index.return=T, decreasing=T)

    #@ Remove start and end points -- hec will not accept connections at start/end of channel
    if(length(station_order$ix)>2){
    #    station_order_inside=station_order$ix[2:(length(station_order$ix)-1)]
    #    weir_pts=weir_pts[station_order_inside,] 

    #    #@ Re-define key variables
    #    st1=as.character(weir_pts$station_name)
    #    st2=gsub('\\*', '', st1) # Remove * symbol
    #    st3=as.numeric(st2) # Now they are numbers

    #    #@ Sort the weir points along the channel
    #    station_order=sort(st3,index.return=T, decreasing=T)
        SINGLE_POINT_FLAG=0 #The weir covers multiple cross-sections
    }else{
        print('##')
        print(paste('WARNING: There is an intersection of', storage_name, ' with the channel.'))
        print('However, it does not contain more than 2 bank points')
        print('I will try to add a single-point lateral weir instead')
        print(' ')
        SINGLE_POINT_FLAG=1 # The weir covers at most 1 cross-section.
        #return(hec_lines)
    }

    #browser()
    tmp_coord=as.character(coordinates(weir_pts)[1,])
    tmp_name=weir_pts$reach_name[1]
    tmp_stat=weir_pts$station_name[1]
    print(' ')
    print('##############################')
    print(paste('Connecting ', storage_name, ' to the river', tmp_name, 'at station', tmp_stat, ' near:'))
    print(tmp_coord)
    print('##############################')
    print(' ') 
    #@ Get coordinates of the weir points, and downstream distances, and make a
    #@ line along them with 3 times as many points. We will get weir elevations
    #@ along this line
    weir_line = coordinates(weir_pts)
    ll=length(weir_line[,1])
    if(ll<2){
        print(' ')
        print('#####################################################')
        print('WARNING: Storage area intersects channel at only one point') 
        print('The relevant point is:')
        print(weir_pts)
        print('No lateral weir will be made here')
        print('#####################################################')
        print(' ')
        return(hec_lines)
    }
    #browser()
    interpolated_length=min(ll+2,80)
    weir_line=cbind(weir_line, c(0, cumsum(weir_pts$downstream_distance[1:(ll-1)])), weir_pts$bank_elev) # Append distance and bank elevation
    # Slightly shorten the downstream distances of weir line: This is a hack to
    # prevent slight overshoots of the weir line beyond the downstream cross-section, which occur for other reasons.
    if(max(weir_line[,3])>0.){
        weir_line[,3] = weir_line[,3]/(max(weir_line[,3]))*(max(weir_line[,3])-3.0)
        if(max(weir_line[,3])<0.0) stop('ERROR: Weir lines are now too short')
    }
    weir_line_xint = approx(weir_line[,3], weir_line[,1], n=interpolated_length) # Interpolate x's
    weir_line_yint = approx(weir_line[,3], weir_line[,2], n=interpolated_length) # Interpolate y's
    weir_line_bank_elev_limit=approx(weir_line[,3], weir_line[,4], n=interpolated_length) # Interpolate the bank elevation
    #@ weir_line = downstream_distance, x_coordinate, y_coordinate
    weir_line=cbind(weir_line_xint$x, weir_line_xint$y, weir_line_yint$y, weir_line_bank_elev_limit$y)

    #@ Sample the raster elevations at points on weir_line
    transect_inds = cbind(rowFromY(lidar_DEM, weir_line[,3]), colFromX(lidar_DEM, weir_line[,2]) )
    weir_elev=lidar_DEM[transect_inds]
    weir_elev=weir_elev + vertical_datum_offset
    #print("############")
    #print(weir_elev)
    if(limit_weir_elevation_by_channel_bank_elevation){
        # Make sure that the weir elevation is >= the channel bank elevation
        weir_elev=pmax(weir_elev, weir_line[,4])
    }
    #print(weir_elev)

    # Ensure that min weir elev is not below the lower_limit_on_lateral_weir_elevations
    if(lower_limit_on_lateral_weir_elevations!=-Inf){
        print(paste('CAREFUL: Forcing the lateral weir elevations to be > lower_limit_on_lateral_weir_elevations ( ', 
                     lower_limit_on_lateral_weir_elevations, ')'))
        weir_elev=pmax(weir_elev, lower_limit_on_lateral_weir_elevations)
    } 
        
    weir_relation=cbind(weir_line[,1], weir_elev)
   
    #@
    #@ Develop text information for write out modified hec_lines
    #@ 

    #@ Define weir station name, just downstream of the most upstream bank station
    #@ Use a different distance for left and right bank points -- trick to
    #@ avoid weirs on both banks with the same station name
    lateral_weir_distance_name=(0.25 + 0.5*(bank2=='R'))*(station_order$x[1]-station_order$x[2])
    #lateral_weir_distance=(0.25 + 0.5*(bank2=='R')) + min(1.0,(station_order$x[1]-station_order$x[2]))
    weir_station_name=station_order$x[1]- lateral_weir_distance_name

    output_text=c()
    line1= paste("Type RM Length L Ch R = 6 ,", round(weir_station_name,3), ",,, ", sep="")
    output_text=c(output_text,line1)
    output_text=c(output_text, 'BEGIN DESCRIPTION:')
    output_text=c(output_text, paste(storage_name, '_to_chan', sep=""))
    output_text=c(output_text, 'END DESCRIPTION:')

    #@ Add timestamp
    edit_time=format(Sys.time(), "%b/%d/%Y %H:%M:%S") # Timestamp like hec-ras
    output_text=c(output_text, paste("Node Last Edited Time=", edit_time,sep=""))

    bankflag=0 + 3*(bank2=='R') # 0 for left bank, 3 for right bank
    output_text=c(output_text, paste("Lateral Weir Pos= ",bankflag, sep=""))


    output_text=c(output_text,
                  paste("Lateral Weir End=                ,                ,        ,",storage_name,sep=""))

    #output_text=c(output_text, paste('Lateral Weir Distance=', round(lateral_weir_distance,3),sep=""))
    output_text=c(output_text, paste('Lateral Weir Distance=', round(1.000,3),sep="")) # Make 1m downstream of upstream section

    if(SINGLE_POINT_FLAG==0){
        lw_tw_text="Lateral Weir TW Multiple XS=-1"
    }else{
        lw_tw_text="Lateral Weir TW Multiple XS=0"
    }
    next_lines=c(lw_tw_text,"Lateral Weir WD=20","Lateral Weir Coef=1.1",
                 "Lateral Weir WSCriteria=-1","Lateral Weir Flap Gates= 0","Lateral Weir Hagers EQN= 0 ,,,,,",
                 "Lateral Weir SS=0.05,0.05,","Lateral Weir Type= 0","Lateral Weir Connection Pos and Dist= 0 ,")

    output_text=c(output_text,next_lines)
    output_text=c(output_text, paste("Lateral Weir SE= ",interpolated_length,sep=""))

    #@ Add weir distance-elevation information
    weir_x=pad_string(as.character(round(weir_relation[,1],3)), 8, pad=" ", justify='right')
    weir_y=pad_string(as.character(round(weir_relation[,2],3)), 8, pad=" ", justify='right')
    weir_text=paste(weir_x,weir_y,sep="")
    weir_text=format_in_rows(weir_text,5)
    output_text=c(output_text, weir_text)

    upstream_station=pad_string(as.character(station_order$x[1]), 8, pad=" ", justify='left')

    output_text=c(output_text, 
                  paste("Lateral Weir HW RS Station=", 
                         upstream_station,",", 
                         -round(lateral_weir_distance_name,3),sep=""))

    output_text=c(output_text, "Lateral Weir TW RS Station=,0")
    output_text=c(output_text, "LW Div RC= 0 ,False,")

    #@
    #@
    #@ Now find the line in 'hec_lines' where we need to insert this weir
    #@ Note that we will be altering hec_lines with every call to this function
    #@ So we need to recompute the locations
    #@


    #@ Identify the start of the channel of interest
    chan_ind=grep('River Reach=', hec_lines) # All channels
    my_chan_ind=grep(weir_pts$reach_name[1], hec_lines[chan_ind]) # Should only match 1

    #@ Compute lower bound of 'lines of interest'
    my_lower_ind=chan_ind[my_chan_ind]
    #@ Compute upper bound of 'lines of interest'
    ll=length(hec_lines)
    if(my_chan_ind<length(chan_ind)){
        my_upper_ind=chan_ind[my_chan_ind+1]
    }else{
        my_upper_ind=ll
    }
  
    #browser() 
    #@ Find the line at the start of the upstream station
    line_pattern=paste('Type RM Length L Ch R = 1 ,',weir_pts$station_name[1],sep="")
    #@ NOTE: Match can contain a * -- need to use 'fixed=TRUE' to get this
    upstream_station_ind=grep(line_pattern, hec_lines[my_lower_ind:my_upper_ind],fixed=TRUE) + my_lower_ind-1

    #@ Find next blank line (search next 300 lines). Identified because it has
    #@ no letters or numbers
    upper_search_ind=min(upstream_station_ind+300, ll)
    blank_loc=grep("[a-z A-z 0-9]", hec_lines[upstream_station_ind:upper_search_ind], invert=TRUE)[1] + upstream_station_ind-1
    #print(blank_loc)
    hec_linestmp=c(hec_lines[1:blank_loc], output_text, " ", hec_lines[(blank_loc+1):ll])

    hec_linestmp
}


################################################################################

usdistfun<-function(point_coords, line_coords, roundoff_tol=1.0e-03){
    # Compute the upstream distance of point='point_coords'=c(x0,y0), along a line
    # defined by 'line_coords' = matrix with 2 columns
    #
    # Also, compute the index of the segment in line_coords that point_coords
    # intersects with. This was most convenient to calculate here, even though
    # its a bit out of place 
    #
    # point_coords should lie on line_coords to within a small tolerence. This
    # tolerance is is related to roundoff_tol in a complex fashion -- read the
    # source if you need to know!
    #
    # NOTE: line_coords is assumed to be ordered from upstream to downstream.
    # Upstream distance is measured from the most downstream point in
    # line_coords
     
    x0=point_coords[1]
    y0=point_coords[2]

    if(length(line_coords[,1])<2) stop('ERROR: line_coords has < 2 points')
    if(length(point_coords)!=2) stop('ERROR: point_coords does not describe a single point')

    # Compute upstream distance along line
    line_coords_rev=line_coords[length(line_coords[,1]):1,]
    us_dist=cumsum( (diff(line_coords_rev[,1])**2 + diff(line_coords_rev[,2])**2)**0.5)
    us_dist=c(0, us_dist)
    us_dist=rev(us_dist) # Indices follow those in line_coords

    output=NA
    connecting_segment=NA
    for(k in 1:(length(us_dist)-1)){
        # Check if x,y lies on this segment of us_dist
        # Do this by computing the change in x and y for 1 unit of movement
        # along the segment, and then taking a segment to x0,y0, computing
        # the same, and comparing
        ds = us_dist[k]-us_dist[k+1]  # This will be positive   
        dy=(line_coords[k+1,2]-line_coords[k,2])/ds # delta y
        dx=(line_coords[k+1,1]-line_coords[k,1])/ds # delta y

        ds2 = (( line_coords[k+1,2]-y0)**2 + (line_coords[k+1,1]-x0)**2)**0.5

        # Quick check to see if we are on possibly on this segment
        if( ds2>ds+roundoff_tol){
            # x0, y0 is not on the segment
            next
        }
        # Allow for degenerate case -- avoid division by ds2
        if(ds2<roundoff_tol){
            output=us_dist[k+1]
            connecting_segment=k
            break
        }

        dy2 = (line_coords[k+1,2]-y0)/ds2
        dx2 = (line_coords[k+1,1]-x0)/ds2

        # Test for equity, allowing for floating point error if ds2 is very small
        if((abs(dy2-dy)<roundoff_tol) && (abs(dx2-dx)<roundoff_tol)){
            # FIXME: roundoff_tol here should be different to the one used previously.
            # Might matter in unusual cases
            output=us_dist[k+1] + ds2
            connecting_segment=k
            break
        }
    }
    

    if(is.na(output) | is.na(connecting_segment)){
        print('output / connecting_segment is NA in usdistfun ( useful_functions.R). Going into browser()')
        browser()
        print('XXX')
        print('XXX')
        print('XXX')
    }
    return(c(output, connecting_segment))
}

