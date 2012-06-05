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
    for(i in 1:(length(chan_ind)-1)){
        offset=chan_ind[i]-1
        #@ Only grep the relevant cross-sections
        xsect_start=grep('XS GIS Cut Line=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        #xsect_end=grep('Node Last Edited Time=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        xsect_end=grep('#Sta/Elev=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset - 1
   
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

        
    }
        #@ Now coerce output to SpatialPolygonsDataFrame
        chan=SpatialPolygons(chan_polygons, proj4string=spatial_proj)

        chan2=SpatialPolygonsDataFrame(chan, data=data.frame(id=seq(1,length(chan_polygons))), match.ID=FALSE)

        chan2
}

#################################################################################################

compute_stage_vol_relation<-function(my_poly,lidar_DEM, vertical_datum_offset, upper_bound_stage=100){
    #@
    #@ Function to compute the stage-volume relation for a region defined by the polygon my_poly
    #@ , on a given lidar_Dem raster  
    #@ Output: A Stage-Volume relation, with the stage values offset by + vertical_datum_offset
    #@

    cell_area=xres(lidar_DEM)*yres(lidar_DEM)/1000. # Area of a pixel -- units should match the units of the output data.
    # In hec-ras, the storage area units are in 1000s of m^3 -- so we need 10
    # lots of 10*10 squares filled in a metre of water to account for one unit.
    # Hence the strange use of units above

    elev_pts=lidar_DEM[my_poly] # Elevation points inside the polygon
    elev_hist=hist(elev_pts[[1]],n=60)

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

    #@ Compute coordinates inside poly
    storage_area_central_coords=coordinates(storage_area)

    #@ Compute bounding coordinates -- got to love the notation
    storage_area_bounding_coords=coordinates(gBoundary(storage_area))[1][[1]][[1]]

    #@ Remove final bounding coordinate (which is a repeat of the first point, and does not feature in hec)
    l=length(storage_area_bounding_coords[,1])
    storage_area_bounding_coords=storage_area_bounding_coords[-l,]

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
    elev_coord=format(as.character(round(elev_vol[,1],2)), width=8,justify='right',trim=T)
    #vol_coord=elev_vol[,2]
    vol_coord=format(as.character(signif(elev_vol[,2],5)),width=8,justify='right',trim=T)

    next_lines=paste(elev_coord,vol_coord,sep="")
    # Get the format right
    next_lines = format_in_rows(next_lines,5)

    output_text=c(output_text,next_lines)

    output_text
    #print(output_text)
}

####################################################################################################################

make_storage_connection_text<-function(store1, store2, name1, name2, lidar_DEM){
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
        
    intersection=gIntersection(store1, store2) # Polygon containing the intersection of the 2 storage areas
    
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
    elev_relation=lidar_DEM[intersection][[1]]
    l=length(elev_relation)
    elev_relation2=sort(elev_relation)
    #@ How to compute the weir length??
    #@ Idea:
    #@ gLength(intersection) = poly boundary length
    #@ If we assume poly is long and thin, then length~ = boundary length /2 
    #@ We could make it of slightly shorter length to be conservative (e.g. boundary_length *1/3 or *5/12)
    weir_length = gLength(intersection)*1/3

    weir_x_vals=seq(0,weir_length, len=l) # X values at which we will get weir elevation points
    l2=min(l,10) # Number of points on the weir in hec-ras

    weir_relation=approx(weir_x_vals, elev_relation2, n=l2) # Weir x - elev relation

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
    for(i in 1:(length(chan_ind)-1)){

        offset=chan_ind[i]-1
        #@ Only grep the relevant cross-sections
        xsect_start=grep('XS GIS Cut Line=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        #xsect_end=grep('Node Last Edited Time=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
        xsect_end=grep('#Sta/Elev=', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset - 1
       
        # Useful to have the 'station' lines as well. Initially I tried to use
        # this as the definition of xsect_start, however, I decided that was
        # not the best. However, we need this line to get the station number
        station_start=grep('Type RM Length L Ch R =', hec_lines[chan_ind[i]: chan_ind[i+1]]) + offset
   
        #@ NOW EXTRACT THE START AND END POINTS OF EACH CROSS-SECTION 

        #@ Predefine variables to store start and end points of
        #@ cross-sections
        coords_start=c()
        coords_end=c()
        reach_name=c()
        station_name=c()
        left_bank_downstream_dist=c()
        right_bank_downstream_dist=c()
        for(j in 1:length(xsect_start)){  # Loop over all xsections
            cutline_text=hec_lines[(xsect_start[j]+1):(xsect_end[j]-1)]
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

        }
       
        output_coords=rbind(output_coords, coords_start)
        output_data=rbind(output_data, 
                          cbind(reach_name, station_name, 
                                rep('L', length(reach_name)), 
                                left_bank_downstream_dist) )
        output_coords=rbind(output_coords, coords_end)
        output_data=rbind(output_data, cbind(reach_name, station_name, 
                                             rep('R', length(reach_name)), 
                                             right_bank_downstream_dist))
         
        #print(cbind(station_name, reach_name))

    }

        # Coerce to spatial points
        output_pts=SpatialPointsDataFrame(coords=output_coords[,1:2], 
                      data=data.frame(reach_name=output_data[,1], station_name=output_data[,2],
                                      bank=output_data[,3], downstream_distance=as.numeric(output_data[,4])),
                      match.ID=FALSE,
                      proj4string=CRS(spatial_proj))
        output_pts
}

##############################################################

make_lateral_weir_text<-function(storage_poly, storage_name, chan_poly, chan_pts, lidar_DEM, hec_lines){
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
  
    #@ Get channel boundary points 
    weir_pts=chan_pts[intersection,]

    #@ Check to ensure we intersect the channel on only 1 bank
    banks=weir_pts$bank
    bank2=union(banks,banks)
    if(length(bank2)>1){
        print(paste('ERROR: Storage polygon,', storage_name, "intersects on both the left and right banks"))
        print(weir_pts)
        stop('FIX THIS BEFORE CONTINUING')
    }

    #@ Sort the bank stations
    st1=as.character(weir_pts$station_name)
    st2=gsub('\\*', '', st1) # Remove * symbol
    st3=as.numeric(st2) # Now they are numbers

    #@ Sort the weir points along the channel
    station_order=sort(st3,index.return=T, decreasing=T)
    weir_pts=weir_pts[station_order$ix,] 

    #@ Get coordinates of the weir points, and downstream distances, and make a
    #@ line along them with 3 times as many points
    weir_line = coordinates(weir_pts)
    ll=length(weir_line[,1])
    interpolated_length=3*ll
    weir_line=cbind(weir_line, c(0, cumsum(weir_pts$downstream_distance[1:(ll-1)]))) # Append distance
    weir_line_xint = approx(weir_line[,3], weir_line[,1], n=interpolated_length) # Interpolate x's
    weir_line_yint = approx(weir_line[,3], weir_line[,2], n=interpolated_length) # Interpolate y's
    #@ weir_line = downstream_distance, x_coordinate, y_coordinate
    weir_line=cbind(weir_line_xint$x, weir_line_xint$y, weir_line_yint$y)

    #@ Sample the raster elevations at points on weir_line
    transect_inds = cbind(rowFromY(lidar_DEM, weir_line[,3]), colFromX(lidar_DEM, weir_line[,2]) )
    weir_elev=lidar_DEM[transect_inds]

    weir_relation=cbind(weir_line[,1], weir_elev)
   
   
    #@
    #@ Develop text information for write out
    #@ 

    #@ Define weir station name, just downstream of the most upstream bank station
    lateral_weir_distance=0.02*(station_order$x[1]-station_order$x[2])
    weir_station_name=station_order$x[1]- lateral_weir_distance

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

    output_text=c(output_text, paste('Lateral Weir Distance=', round(lateral_weir_distance,3),sep=""))

    next_lines=c("Lateral Weir TW Multiple XS=-1","Lateral Weir WD=20","Lateral Weir Coef=1.1",
                 "Lateral Weir WSCriteria=-1","Lateral Weir Flap Gates= 0","Lateral Weir Hagers EQN= 0 ,,,,,",
                 "Lateral Weir SS=0.05,0.05,","Lateral Weir Type= 0","Lateral Weir Connection Pos and Dist= 0 ,")

    output_text=c(output_text,next_lines)
    output_text=c(output_text, paste("Lateral Weir SE= ",interpolated_length,sep=""))

    #@ Add weir distance-elevation information
    weir_x=pad_string(as.character(signif(weir_relation[,1],7)), 8, pad=" ", justify='right')
    weir_y=pad_string(as.character(signif(weir_relation[,2],7)), 8, pad=" ", justify='right')
    weir_text=paste(weir_x,weir_y,sep="")
    weir_text=format_in_rows(weir_text,5)
    output_text=c(output_text, weir_text)

    upstream_station=pad_string(as.character(station_order$x[1]), 8, pad=" ", justify='left')

    output_text=c(output_text, 
                  paste("Lateral Weir HW RS Station=", 
                         upstream_station,",", 
                         -round(lateral_weir_distance,3),sep=""))

    output_text=c(output_text, "Lateral Weir TW RS Station=,0")
    output_text=c(output_text, "LW Div RC= 0 ,False,")

    #@
    #@
    #@ Now find the line in 'hec_lines' where we need to insert this weir
    #@
    #@

    chan_ind=grep('River Reach=', hec_lines)

}
