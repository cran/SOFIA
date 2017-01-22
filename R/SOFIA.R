SOFIA<-function(data,chromoConfiguration=NULL,dataColorFlag=FALSE,dataColor=NULL,plotType=NULL,plotColor=NULL,markerSize=NULL,plotLocation=NULL,plotBackground=NULL,plotImportance=NULL,plotOrientation=NULL,density=NULL,linksFlag=FALSE,linkColor='red',linkGeometry=c(.2,.2),linkRadius=c(.7,.7),blocksFlag=FALSE,blocksData=NULL,tilesFlag=FALSE,tilesData=NULL,tilesLocation=NULL,blocksColor=rbind(c('blue','red'),c('green','purple')),blocksLocation=c(.5,.01),gaps=NULL,ideogramThickness=20,generalPlotConfFlag=FALSE,generalPlotConf=NULL,chrPrefixFont='upper',chrPrefix='LG',tickSeparation=50,tickSuffix='cM',circosLocation=NULL,returnConf=FALSE,circosDisplay=FALSE,figureDisplay=TRUE,deleteData=TRUE,runCircos=TRUE,confName='circos.conf'){
  
  
  namename<-confName
  
  if (is.null(circosLocation)){
    cat('\ncircosLocation was not defined, please provide the directory where /bin/circos is contained. This argument is requiered.\n')
  }else{
  kr<-list.files()
  kr2<-match(c('circos.png','circos.svg'),kr)
  
  for (kt in 1:2){
    if (is.na(kr2[kt])==FALSE){
      cmd <- paste0('rm ', getwd(),'/',kr[kr2[kt]])
      system(cmd)
    }
  }
  
  
  if(runCircos){
  cat(paste0('\nIf no errors are produced, you can find circos.png and circos.svg (ediatable version) in ',getwd(),'\n'))
  }else{
  cat('\nSince runCircos=FALSE, Circos was not invoked. Your configuration and data files can be find in /circos/bin/ and circos/data/. If you want SOFIA to run Circos, set this parameter TRUE (or just remove it)')  
  }
  
	nPheno<-ncol(data)-4
	if (nPheno==0 & linksFlag==FALSE){
		cat('\nYou did not provide data to plot. If you only want to plot links between two or more maps, please set linksFlag=TRUE and populate the corresponing parameters\n')
	}
	
	randomColors<-c('red','blue','green','purple','orange','yellow','dred','dblue','dgreen','dpurple','dorange','dyellow','lred','lblue','lgreen','lpurple','lorange','lyellow','vdred','vdblue','vdgreen','vdpurple','vdorange','vdyellow','vlred','vlblue','vlgreen','vlpurple','vlorange','vlyellow')
	
	if (is.null(plotType) & nPheno>0){
		plotType<-numeric()
		for (i in 1:nPheno){
			if (is.numeric(data[,4+i])){
			plotType<-c(plotType,'scatter')
		}else{
			plotType<-c(plotType,'text')
		}
		}
	}
	
	if (is.null(plotColor) & nPheno>0){
		plotColor<-numeric()
		for (i in 1:nPheno){
			plotColor<-c(plotColor,randomColors[sample(1:15,1)])
		
			
		}
		}
	
	
	
	if (is.null(plotLocation) & nPheno>0){
		stPosition<-.99
		sizeF<-(.8/nPheno)
		plotLocation<-data.frame(r0=rep(NA,nPheno),r1=rep(NA,nPheno))
		for (i in 1:nPheno){
			plotLocation$r0[i]<-stPosition-sizeF
			plotLocation$r1[i]<-stPosition
			stPosition<-stPosition-sizeF
		}
	}
	
	
	
	if (blocksFlag & is.null(blocksData)==TRUE){
		stop('\n You set blocksFlag=TRUE but no data was provided in argument blocksData\n')
	}
	
	if (tilesFlag & is.null(tilesData)==TRUE){
		stop('\n You set tileFlag=TRUE but no data was provided in argument tilesData\n')
	}
	
	
	if (is.null(plotImportance)){
		plotImportance<-rep(1,nPheno)
	}
	
	if (is.null(plotBackground)){
		plotBackground<-data.frame(backgroundShow=rep(FALSE,nPheno),axisShow=rep(FALSE,nPheno))
	}
	
	
	
	
	if (is.null(plotOrientation)){
		plotOrientation<-rep('out',nPheno)
	}
	
	if (is.null(markerSize)){
		markerSize<-rep(10,nPheno)
	}
	
	
	
	# universal variables
	options(scipen=999)
	multiplicador<-1000


	chrMapper<-as.character(paste0(data$map,'x',data$chr))
	numberMaps<-length(unique(data$map))
	numberChromosomes<-length(unique(chrMapper))
	numberChromosomesMap<-numeric()
	for (i in 1:numberMaps){
		f<-which(data$map==i)
		numberChromosomesMap<-c(numberChromosomesMap,length(unique(data$chr[f])))
	}
	
	
	# making a chromoMap
	superMap<-unique(chrMapper)
	superMap<-rbind(superMap,1:numberChromosomes)
	chrMapperReal<-chrMapper
	for (i in 1:numberChromosomes){
		f<-which(chrMapper==superMap[1,i])
		chrMapperReal[f]<-rep(i,length(f))
	}
	
	
	if (tilesFlag){
	chrMapper_t<-as.character(paste0(tilesData$map,'x',tilesData$chr))
	chrMapperReal_t<-chrMapper_t
	for (i in 1:numberChromosomes){
		f<-which(chrMapper_t==superMap[1,i])
		if (length(f)>0){
		chrMapperReal_t[f]<-rep(i,length(f))
		
	}
	}
	}
	

	if (substring(circosLocation,nchar(circosLocation))!='/'){
		circosLocation<-paste0(circosLocation,'/')
	}
	
	
	w1<-list.files(paste0(circosLocation,"bin/"),pattern='circos')
	if (length(w1[w1=='circos'])){
		cat(paste0('\ncircos was found in ',paste0(circosLocation,'bin/'),'\n\n','starting SOFIA...\n\n\n'))
	}else{
		cat(paste0('\nERROR: /bin/circos was not found in the provided directory (',paste0(circosLocation,"bin/"),')\nPlease try again...\n\n\nbye!!!\n\n\n'))
	}
	
	
	
	## making ticks for ideograms	
	tiks<-c('show_ticks=yes','show_tick_labels=yes','show_grid=no','<ticks>','radius=dims(ideogram,radius_outer)','label_offset=5p','label_size=8p','multiplier=1e-3','color=black','<tick>',paste0('spacing=',tickSeparation,'u'),'size=8','thickness=2p','show_label=yes','label_size=20',paste0('suffix=',tickSuffix),paste0('grid_start=.9r'),'grid_end=1r+45p','grid_color=vdgrey','grid_thickness=1p','grid=yes','</tick>','</ticks>')
	write(tiks, file = paste0(circosLocation,'/bin/tiks.conf'))
	 
	 
	 
	# making ideograms
	chromoSizes<-numeric()
	
	for (i in 1:numberChromosomes){
		f<-which(chrMapper==unique(chrMapper)[i])
		chromoSizes<-c(chromoSizes,max(data$pos[f]))
	}
	 
	 
	 
	
	
	# making defaults for chromoConfiguration
	if (is.null(chromoConfiguration)){
		chromoConfiguration<-data.frame()
	}
	if (is.null(chromoConfiguration$color)){
	ideoColors<-numeric()
	for (i in 1:numberMaps){
		ideoColors<-c(ideoColors,rep(randomColors[i],numberChromosomesMap[i]))
	}
}else{
		ideoColors<-chromoConfiguration$color
	}
	
	
	if (is.null(chromoConfiguration$order)){
	xOrder<-1:numberChromosomes
}else{
	zz<-paste0(chromoConfiguration$map,'x',chromoConfiguration$order)
	zz2<-match(zz,superMap[1,])
		xOrder<-superMap[2,zz2]
	}
	
	if (is.null(chromoConfiguration$rev)){
	xRev<-rep(FALSE,numberChromosomes)
}else{
		xRev<-chromoConfiguration$rev
	}
	
	if (is.null(chromoConfiguration$radius)){
	xRad<-rep(1,numberChromosomes)
}else{
		xRad<-chromoConfiguration$radius
	}
	
	 
	# making ideogram file 
	ideo1<-numeric()
	for (i in 1:numberChromosomes){
			ideo1<-c(ideo1,paste(paste('chr - c',i,sep=''),superMap[1,i],0,round(chromoSizes[i]*multiplicador)+20,ideoColors[i],sep=' '))
	}

	write(ideo1, file = paste0(circosLocation,"data/ideo1.txt"))
	
	
	# starting confFile
	confFile<-c('<<include etc/colors_fonts_patterns.conf>>','<<include bin/tiks.conf>>','karyotype = data/ideo1.txt',paste0('chromosomes_units=',multiplicador),'chromosomes_display_default=yes',
	'<<include ideogramLalo.conf>>','<image>','radius* = 1500p','<<include etc/image.conf>>','</image>',paste0('chromosomes_radius=',paste0(paste0('c',1:numberChromosomes,':',xRad,'r'),collapse=",")),paste0('chromosomes_order=',paste0(paste0('c',xOrder),collapse=",")),if(length(which(xRev==TRUE))>0){paste0('chromosomes_reverse=',paste0('c',which(xRev==TRUE),collapse=","))},'<<include etc/housekeeping.conf>>')


	# making ideogram conf file
	
	spacesIdeo<-numeric()
	if (is.null(gaps)==FALSE){
		mybreaks<-cbind(paste0(gaps$mapA,'x',gaps$chrA),paste0(gaps$mapB,'x',gaps$chrB))
		for (i in 1:nrow(gaps)){
			myb<-c(which(mybreaks[i,1]==superMap[1,]),which(mybreaks[i,2]==superMap[1,]))
			spacesIdeo<-c(spacesIdeo,paste0('<pairwise c',myb[1],' c',myb[2],'>'),paste0('spacing=',gaps$length[i],'u'),'</pairwise>')
		}
	}
	
	ideoFile<-c('<ideogram>','show=yes','<spacing>','default=8u',if(length(spacesIdeo)==0){''}else{spacesIdeo},'</spacing>','label_color=black','radius=.88r',paste0('thickness=',ideogramThickness,'p'),'fill = yes','show_label=yes','label_radius=1.08r','label_size=20','label_parallel=yes',paste0('label_case=',chrPrefixFont),'show_bands=yes','fill_bands=yes','</ideogram>')

	
	write(ideoFile, file = paste0(circosLocation,"bin/ideogramLalo.conf"))


	if (linksFlag){
		if (numberMaps==1){stop('\nYou can not make links with only one map\n')}
				combos<-combn(1:numberMaps,2)
				link1<-numeric()
				cont<-1
			
			for (g in 1:ncol(combos)){
				matrixX<-numeric()
				arco1<-which(data$map==combos[1,g])
				arco2<-which(data$map==combos[2,g])
				x1<-match(data$locus[arco1],data$locus[arco2])
				for (i in 1:length(x1)){
					if (is.na(x1[i])==F){
						link1[cont]<-paste(paste0('c',chrMapperReal[arco1[i]]),round(multiplicador*data$pos[arco1[i]]),round(multiplicador*data$pos[arco1[i]]+20),paste0('c',chrMapperReal[arco2[x1[i]]]),round(multiplicador*data$pos[arco2[x1[i]]]),round(multiplicador*data$pos[arco2[x1[i]]]+20),if(linkColor=='chr'){paste0('color=',ideoColors[as.numeric(chrMapperReal[arco1[i]])])}else{paste0('color=',linkColor)})
						cont<-cont+1
					}
				}
			}

				write(link1, file = paste0(circosLocation,"data/x_link1.txt"))


				# adding this to configuration file
	
				confFile<-c(confFile,'<links>',paste0('radius1=',linkRadius[1],'r'),paste0('radius2=',linkRadius[2],'r'),paste0('bezier_radius=',linkGeometry[1],'r'),paste0('crest=',linkGeometry[2]),'<link segdup>','show=yes','file=data/x_link1.txt','record_limit=30000','</link>','</links>')				
		}else{
			cat('\nlinksFlag is FALSE. If you have at least two maps with at least one locus in common you can add links to your plot. Please set linksFlag=TRUE and populate the corresponding parameters\n')
		}
		
		
		
	dedotes<-FALSE
	
	
	if (blocksFlag){
	contx2<-1
	badIdea<-rep(c(1,2),1000)
	colorLD<-blocksColor
	block1<-numeric()
	blockData<-list()

	for (j in 1:ncol(blocksData)){  #number of invidivduals as columns in data3. Rows must be the markers (=# markers in map)
		contX<-1
		ldBlock<-numeric()
		for (i in 1:numberMaps){
			vat1<-which(data$map==unique(data$map)[i])
			
			
			
			for (k in 1:length(unique(data$chr[vat1]))){
				
					f<-which(data$chr[vat1]==unique(data$chr[vat1])[k])
					f0<-c(as.numeric(blocksData[vat1[f],j]),3)
					
					
					
				if (f0[1]==0){
					cont<-1
				}else{
					cont<-2
				}
				cont2<-1
				block1<-numeric()
				for (m in 2:length(f0)){
					if (f0[m]==f0[m-1]){
						cont2<-cont2
					}else{
						block1<-rbind(block1,c(cont2,m,badIdea[cont]))
						cont<-cont+1
						#cont0<-cont0+1
						cont2<-m
					}

				}
				block1[nrow(block1),2]<-block1[nrow(block1),2]-1
				for (n in 1:nrow(block1)){
						ldBlock[contX]<-paste(paste0('c',unique(chrMapperReal[vat1[f]])),round(multiplicador*data$pos[vat1[f[block1[n,1]]]]),round(multiplicador*data$pos[vat1[f[block1[n,2]]]]),data$pos[vat1[f[block1[n,1]]]],paste0('color=',colorLD[i,block1[n,3]]))
					contX<-contX+1
				}
			}

		}
		blockData[[contx2]]<-ldBlock
		contx2<-contx2+1	
	}


	trStart<-0
	for (i in 1:ncol(blocksData)){
		fileX<-paste0(paste0(circosLocation,"data/tr_auto_"),trStart,".txt")
		write(blockData[[i]], file = fileX)
		trStart<-trStart+1
	}

	## making plotConf

	plotConf<-c('<plot>','<<include file.conf>>','<<include location.conf>>','type=heatmap','</plot>')
	write(plotConf, file = paste0(circosLocation,"bin/plot.conf"))

	## making r0r1.conf

	r1r0Conf<-c('r0=eval(sprintf("%fr",conf(track_start)-counter(plot)*conf(track_step)))','r1=eval(sprintf("%fr",conf(track_start)+conf(track_width)-counter(plot)*conf(track_step)))','orientation = eval( counter(plot) % 2 ? "in" : "out" )')
	write(r1r0Conf, file = paste0(circosLocation,"bin/location.conf"))

	## making file.conf

	fileConf<-'file = data/tr_auto_counter(plot).txt'
	write(fileConf, file = paste0(circosLocation,"bin/file.conf"))

	## updating confFile
	if (generalPlotConfFlag){
	  confFile<-c(confFile,paste0('track_width=',blocksLocation[2]),paste0('track_start=',blocksLocation[1]),paste0('track_step=',blocksLocation[2]+0.0005),'<plots>',generalPlotConf,rep('<<include plot.conf>>',ncol(blocksData)))
	}else{
	  confFile<-c(confFile,paste0('track_width=',blocksLocation[2]),paste0('track_start=',blocksLocation[1]),paste0('track_step=',blocksLocation[2]+0.0005),'<plots>',rep('<<include plot.conf>>',ncol(blocksData)))
	 	}
	dedotes<-TRUE
}

	
	if(dedotes==FALSE){
	  if (generalPlotConfFlag){
	    confFile<-c(confFile,'<plots>',generalPlotConf)
	  }else{
			confFile<-c(confFile,'<plots>')
	  }
		}
		
		if (nPheno>0){
			for (ii in 1:length(plotType)){
				if (plotType[ii]=='scatter'){
					
					f<-which(is.na(data[,4+ii])==FALSE)
					track2<-paste(paste0('c',chrMapperReal[f]),round(multiplicador*data$pos[f]),round(multiplicador*data$pos[f]+1),data[f,4+ii],if(dataColorFlag){paste0('color=',dataColor[f,ii])}else{if(plotColor[ii]=='chr'){paste0('color=',ideoColors[as.numeric(chrMapperReal[f])])}else{''}})
							
						
			
					write(track2, file = paste0(paste0(circosLocation,'data/tr_2_'),ii,'.txt'))
					
				
					
					backAx<-c('<backgrounds>','show=data','<background>',if(plotBackground$backgroundShow[ii]){paste0('color=',plotBackground$backgroundColor[ii])}else{''},'</background>','</backgrounds>','<axes>','<axis>',if(plotBackground$axisShow[ii]){c('color=vdgrey',paste0('spacing=',plotBackground$axisSep[ii]))}else{''},'</axis>','</axes>')
					
					confFile<-c(confFile,'<plot>','type=scatter',paste0('file=data/tr_2_',ii,'.txt'),paste0('z=',plotImportance[ii]),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),paste0('orientation=',plotOrientation[ii]),paste0('glyph=','circle'),paste0('glyph_size=',markerSize[ii]),paste0('color=',plotColor[ii]),paste0('color_log_scale=','1'),backAx,'</plot>')
			
				}
				if (plotType[ii]=='heatmap'){
														
					f<-which(is.na(data[,4+ii])==FALSE)
					track3<-paste(paste0('c',chrMapperReal[f]),round(multiplicador*data$pos[f]),round(multiplicador*data$pos[f]+markerSize[ii]),data[f,4+ii],if(dataColorFlag){paste0('color=',dataColor[f,ii])}else{if(plotColor[ii]=='chr'){paste0('colorX=',ideoColors[as.numeric(chrMapperReal[f])])}else{''}})
							
						
			
					write(track3, file = paste0(paste0(circosLocation,'data/tr_3_'),ii,'.txt'))
					
				
					confFile<-c(confFile,'<plot>','type=heatmap',paste0('file=data/tr_3_',ii,'.txt'),paste0('z=',plotImportance[ii]),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),paste0('color=',plotColor[ii]),paste0('color_log_scale=','1'),'</plot>')
			
				}
				if (plotType[ii]=='heatmap_interval'){
					myCC<-c('white_a5',plotColor[ii])
				contx2<-1
				badIdea<-rep(c(1,2),1000)
				block1<-numeric()
				blockData<-list()

				j<-ii+4
					contX<-1
					ldBlock<-numeric()
					for (i in 1:numberMaps){
						vat1<-which(data$map==unique(data$map)[i])
			
						for (k in 1:length(unique(data$chr[vat1]))){
				
								f<-which(data$chr[vat1]==unique(data$chr[vat1])[k])
								f0<-c(as.numeric(data[vat1[f],j]),3)
					
					
					
							if (f0[1]==0){
								cont<-1
							}else{
								cont<-2
							}
							cont2<-1
							block1<-numeric()
							for (m in 2:length(f0)){
								if (f0[m]==f0[m-1]){
									cont2<-cont2
								}else{
									block1<-rbind(block1,c(cont2,m,badIdea[cont]))
									cont<-cont+1
									cont2<-m
								}

							}
							block1[nrow(block1),2]<-block1[nrow(block1),2]-1
							
							for (n in 1:nrow(block1)){
									ldBlock[contX]<-paste(paste0('c',unique(chrMapperReal[vat1[f]])),round(multiplicador*data$pos[vat1[f[block1[n,1]]]]),round(multiplicador*data$pos[vat1[f[block1[n,2]]]]),data$pos[vat1[f[block1[n,1]]]],paste0('color=',myCC[block1[n,3]]))
								contX<-contX+1
							}
						}

					}
					blockData[[contx2]]<-ldBlock
					contx2<-contx2+1	
				


				
				write(blockData[[1]], file = paste0(circosLocation,'data/tr_40_',ii,'.txt'))
				
				confFile<-c(confFile,'<plot>','type=heatmap',paste0('file=data/tr_40_',ii,'.txt'),paste0('z=',plotImportance[ii]),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),paste0('color_log_scale=','1'),'</plot>')
				
				
			}
				if (plotType[ii]=='line'){
														
					f<-which(is.na(data[,4+ii])==FALSE)		
					track4<-paste(paste0('c',chrMapperReal[f]),round(multiplicador*data$pos[f]),round(multiplicador*data$pos[f]+1),data[f,4+ii],if(dataColorFlag){paste0('fill_color=',dataColor[f,ii])}else{if(plotColor[ii]=='chr'){paste0('fill_color=',ideoColors[as.numeric(chrMapperReal[f])])}else{''}})
							
						
			
					write(track4, file = paste0(paste0(circosLocation,'data/tr_4_'),ii,'.txt'))
					
				
					
					backAx<-c('<backgrounds>','show=data','<background>',if(plotBackground$backgroundShow[ii]){paste0('color=',plotBackground$backgroundColor[ii])}else{''},'</background>','</backgrounds>','<axes>','<axis>',if(plotBackground$axisShow[ii]){c('color=vdgrey',paste0('spacing=',plotBackground$axisSep[ii]))}else{''},'</axis>','</axes>')
					
					confFile<-c(confFile,'<plot>','type=line',paste0('file=data/tr_4_',ii,'.txt'),paste0('z=',plotImportance[ii]),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),paste0('orientation=',plotOrientation[ii]),paste0('thickness=',markerSize[ii]),paste0('fill_color=',plotColor[ii]),paste0('color_log_scale=','1'),backAx,'</plot>')
			
				}
				if (plotType[ii]=='text'){
														
					f<-which(is.na(data[,4+ii])==FALSE)
							
					track5<-paste(paste0('c',chrMapperReal[f]),round(multiplicador*data$pos[f]),round(multiplicador*data$pos[f]+1),data[f,4+ii])
							

					write(track5, file = paste0(paste0(circosLocation,'data/tr_5_'),ii,'.txt'))
										
					confFile<-c(confFile,'<plot>','type=text',paste0('file=data/tr_5_',ii,'.txt'),paste0('z=',plotImportance[ii]),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),paste0('color=',plotColor[ii]),'show_links=yes','link_dims=0p,5p,5p,5p,0p','link_thickness=1p','link_color=black',paste0('label_size=',markerSize[[ii]],'p'),'label_font=default','padding=0p','rpadding=0p','label_snuggle=yes','max_snuggle_distance=5r','snuggle_tolerance=0.25r','snuggle_sampling=2','</plot>')
			
				}
				if (plotType[ii]=='glyphs'){
														
					f<-which(is.na(data[,4+ii])==FALSE)
							
					track6<-paste(paste0('c',chrMapperReal[f]),round(multiplicador*data$pos[f]),round(multiplicador*data$pos[f]+1),rep('O',length(f)),paste0('color=',data[f,4+ii]))
					backAx<-c('<backgrounds>','show=data','<background>',if(plotBackground$backgroundShow[ii]){paste0('color=',plotBackground$backgroundColor[ii])}else{''},'</background>','</backgrounds>')
					
					write(track6, file = paste0(paste0(circosLocation,'data/tr_6_'),ii,'.txt'))
										
					confFile<-c(confFile,'<plot>','type=text',paste0('file=data/tr_6_',ii,'.txt'),paste0('z=',plotImportance[ii]),'label_font=glyph','padding=-0.1r','rpadding=0p',paste0('label_size=',markerSize[ii],'p'),paste0('r0=',plotLocation$r0[ii],'r'),paste0('r1=',plotLocation$r1[ii],'r'),backAx,'</plot>')
			
				}
			}
			
		}
			# more plots
			if (tilesFlag){
																			
				track6<-paste(paste0('c',chrMapperReal_t),round(multiplicador*tilesData$pos1),round(multiplicador*tilesData$pos2),paste0('color=',tilesData$color))
						
				write(track6, file = paste0(paste0(circosLocation,'data/tr_6_'),ii,'.txt'))
									
				confFile<-c(confFile,'<plot>','type=tile',paste0('file=data/tr_6_',ii,'.txt'),'z=10000',paste0('r0=',tilesLocation[1],'r'),paste0('r1=',tilesLocation[2],'r'),'orientation = out','layers= 25','margin=0.02u','thickness=15','padding=8','stroke_thickness=1','</plot>')
		
			}
			
			if (is.null(density)==FALSE){
			if (density$show){
				track7<-numeric()
				for (ch in 1:numberChromosomes){
						f<-which(unique(chrMapperReal)[ch]==chrMapperReal)
						if (length(f)>0){
						h<-hist(data$pos[f],breaks=density$bins,plot=F)
						h1<-h$breaks[-length(h$breaks)]
						h2<-h$breaks[-1]
						h2[length(h2)]<-chromoSizes[ch]
						track7<-c(track7,paste(paste0('c',unique(chrMapperReal)[ch]),round(multiplicador*h1),round(multiplicador*h2),h$counts))
					
				}
				}
				
						
				write(track7, file = paste0(paste0(circosLocation,'data/tr_7_'),ii,'.txt'))
									
				
				backAx<-c('<backgrounds>','<background>',if(density$backgroundShow){paste0('color=',density$backgroundColor)}else{''},'</background>','</backgrounds>','<axes>','<axis>',if(density$axisShow){c('color=vdgrey',paste0('spacing=',density$axisSep))}else{''},'</axis>','</axes>')
				
				confFile<-c(confFile,'<plot>','type=histogram',paste0('file=data/tr_7_',ii,'.txt'),paste0('z=',density$importance),paste0('r0=',density$r0,'r'),paste0('r1=',density$r1,'r'),paste0('orientation=',density$orientation),paste0('thickness=',density$thickness),paste0('fill_color=',density$color),backAx,'</plot>')
				
			}
		}
			
			
		confFile<-c(confFile,'</plots>')
		write(confFile, file = paste0(circosLocation,paste0("bin/",namename)))
	
  if(runCircos){
	
	if (circosDisplay==FALSE){
		arg1 <- paste0(circosLocation,"bin/circos -silent -conf ",namename)
		cmd <- paste("perl", arg1)
		system(cmd)
	}else{
		arg1 <- paste0(circosLocation,"bin/circos -conf ",namename)
		cmd <- paste("perl", arg1)
		system(cmd)
	}
		
		
	if (returnConf){return(confFile)}
	if (deleteData){
		 cmd <- paste0('rm ', circosLocation,'/data/*.txt')
		 system(cmd)
	}
	if (figureDisplay){	
	if (is.na(kr2[1])==FALSE){
	  tra <- png::readPNG(paste0(getwd(),'/circos.png'))
	  grid::grid.raster(tra)
	}
	}
		
  }	
		
  }
}
