# Utility functions for hec-help

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

    library(sp)
    library(rgdal)

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

compute_stage_vol_relation<-function(my_poly,lidar_DEM, vertical_datum_offset){
    #@
    #@ Function to compute the stage-volume relation for a region defined by the polygon my_poly
    #@ , on a given lidar_Dem.  
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
    upper_bound=100. # MAX value at which we want the volume. This is a trick so that hec-ras never exceeds the range
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

    #@ Compute coordinates inside poly
    storage_area_central_coords=coordinates(storage_area)
    #@ Compute bounding coordinates
    storage_area_bounding_coords=coordinates(gBoundary(storage_area))[1][[1]][[1]]
    #@ Remove final bounding coordinate (which is a repeat of the first point, and does not feature in hec)
    l=length(storage_area_bounding_coords[,1])
    storage_area_bounding_coords=storage_area_bounding_coords[-l,]


    output_text=c() # Initialise output_text

    # Name related character strings
    name_nonempty=as.character(name)
    name_fil=rep(" ", 16-nchar(name_nonempty))
    name_fil=paste(name_fil,collapse='')

    # Coordinate character strings
    coord1=format(as.character(signif(storage_area_central_coords[1,1],13)),width=15,justify='none',trim=T)
    coord2=format(as.character(signif(storage_area_central_coords[1,2],13)),width=15,justify='none',trim=T)

    # Build first line of output text
    first_line=paste("Storage Area=", name_nonempty, name_fil, ",", coord1,",", coord2, sep="")
    #print(first_line)
    output_text=c(output_text, first_line)

    # Build the second line of output text
    second_line=paste('Storage Area Surface Line=', l-1)
    output_text=c(output_text,second_line)

    # Build boundary output coords
    coord1=format(as.character(signif(storage_area_bounding_coords[,1],16)),width=17,justify='none',trim=T)
    coord2=format(as.character(signif(storage_area_bounding_coords[,2],16)),width=17,justify='none',trim=T)
  
    # Forcibly pad these - it is critical to get the width right
    n1=nchar(coord1)
    for(i in 1:length(n1)){
        if(n1[i]<16){
            coord1[i]=paste(paste(rep(" ", 16-n1[i]),collapse=""), coord1[i],sep="")
        }
    }
    
    next_lines=paste(coord1,coord2,sep="")
    #for(i in 1:length(next_lines)){
    #    if(nchar(next_lines[i])<32){
    #        
    #    }
    #
    #}
    output_text=c(output_text,next_lines)

    # Add storage area type info
    output_text=c(output_text, c("Storage Area Type= 1", "Storage Area Area=","Storage Area Min Elev="))

    # Add volume - elev information
    l=length(elev_vol[,1])
    output_text=c(output_text,paste('Storage Area Vol Elev=', l) )
    elev_coord=format(as.character(round(elev_vol[,1],2)), width=8,justify='right',trim=T)
    #vol_coord=elev_vol[,2]
    vol_coord=format(as.character(signif(elev_vol[,2],5)),width=8,justify='right',trim=T)

    next_lines=paste(elev_coord,vol_coord,sep="")

    l=length(next_lines)
    l1= (5-l%%5)%%5 # Number of lines we need to add to make the number of next_lines elements a multiple of 5
    next_lines=c(next_lines, rep(" ", l1))

    next_lines2=matrix(next_lines,ncol=5,byrow=T)
    
    next_lines3=paste(next_lines2[,1], next_lines2[,2], next_lines2[,3], next_lines2[,4], next_lines2[,5], sep="")

    output_text=c(output_text,next_lines3)

    print(output_text)
}
