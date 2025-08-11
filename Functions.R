suppressPackageStartupMessages({
  library(shiny)
  library(SummarizedExperiment)
  library(ggplot2)
  library(dplyr)
  library(Chicago)
  library(cowplot)
  library(stringr)
  library(igraph)
  library(RColorBrewer)
  library(rtracklayer)
  library(netdiffuseR)
  library(IRanges)
  library(tibble)
  library(tidyr)
})

# --- Load Data (CHANGE THESE PATHS) ---
# For deployment, place these files in a 'data' folder and use relative paths
# e.g., file.path("data", "TadBorders.txt")
TadBorders <- read.delim("data/01_TADBoundaris_CNCCE_10.5_StrongBoundaries_250kb.txt")
Prom <- readRDS("data/01_PromoterCategories_FnpMxMdPa2_E10.5_E11.5_genotype_Progenitors_GaussianModel_Test4.rds")
Enh <- readRDS("data/01_EnhancerCategories_FnpMxMdPa2_E10.5_E11.5_genotype_Progenitors_GaussianModel_Test4_after.rds")

# Assuming 'SE', 'E11.5', and 'E10.5' are SummarizedExperiment objects
# You'll need to load these as well.
SE <- readRDS("path/to/your/SE.rds")
E11.5 <- readRDS("path/to/your/E11.5.rds")
E10.5 <- readRDS("path/to/your/E10.5.rds")

########
#Transparency function
#############


Alpha=function(col,fct){
  coul=vector()
  for (i in 1:length(col)) {
    coul[i]=alpha(col[i],alpha=fct[i])
  }
  return(coul)
}


########
#Network Tables
#############

.NetworkTable= function(data=SE,single=FALSE,gene="Hoxa2",Celltype="FNP",Genotype="WT",Time="E11.5",P_P=TRUE,ALL=FALSE,AceGradient=FALSE) {
  
  rdat <- rowData(data)
  score <- assay (data)
  Experiment <- colData(data)
  
  if(Genotype=="WT") {
    CellTypes=subset(Experiment,celltype==Celltype & genotype==tolower(Genotype) & time==Time, select="sample")
  }else{CellTypes=subset(Experiment,celltype==Celltype & genotype==Genotype & time==Time, select="sample")}
  if (single ) { 
    if(P_P) { 
      data.bait <- data[rdat$from.baitSymbol==gene,]
      data.oe <- data[!is.na(rdat$to.oePromoterSymbol),]
      data.oe <- data.oe[rowData(data.oe)$to.oePromoterSymbol==gene,]
      data <- rbind(data.bait,data.oe)
      rdat <- rowData(data)
      score <- assay(data)
      Experiment <- colData(data) 
      
    }else{
      
      data <- data[rdat$from.baitSymbol==gene,] 
      rdat <- rowData(data)
      score <- assay(data)
      Experiment <- colData(data)          
    }
  } 
  
  rdat$chr=rdat$from.X.seqnames
  rdat$start=rowMins(as.matrix(rdat[,c("from.X.start","to.X.start")]),useNames=FALSE)
  rdat$end=rowMaxs(as.matrix(rdat[,c("from.X.end","to.X.end")]),useNames = FALSE)
  Overlap=GRanges(as.data.frame(rdat)%>%select(chr,start,end)) %over% GRanges(TadBorders)  
  rdat$tad="Intra"
  if(sum(Overlap)!=0){rdat[Overlap,]$tad="Inter"}
  
  if (AceGradient){ Col=paste0("H3K27ac_",Celltype,"_",Time,"_",Genotype)
  State=paste0("State_",Celltype,"_",Time,"_",Genotype)
  Prom_df<- as.data.frame(assay(Prom[rowData(Prom)[,State]=="active" & !is.na(rowData(Prom)[,State]) ,]))%>%
    dplyr::select(!!Col)%>%
    dplyr::mutate( !!as.symbol(Col) := replace(!!as.symbol(Col), !!as.symbol(Col) < 0, 0))%>%
    dplyr::mutate(quantile_H3K27ac:= ntile(., 5))
  
  Enh_df<-as.data.frame(assay(Enh[rowData(Enh)$class=="non-promoter" & !is.na(rowData(Enh)$class) &
                                    rowData(Enh)[,State]=="active" & !is.na(rowData(Enh)[,State]) ,]))%>%
    dplyr::select(!!Col)%>%
    dplyr::mutate( !!as.symbol(Col) := replace(!!as.symbol(Col), !!as.symbol(Col) < 0, 0))%>%
    dplyr::mutate(quantile_H3K27ac:= ntile(., 5))
  }	
  
  
  if (P_P) { #PP Networks
    
    Sample <- cbind(rdat,score) %>%
      as.data.frame() %>%
      dplyr::filter(type.oe=="promoter")%>%
      dplyr::select("from.baitSymbol","to.oePromoterSymbol",
                    starts_with("category.bait.") & contains(Time) & contains(Celltype) & contains(Genotype),
                    starts_with("category.oeProm.") & contains(Time) & contains(Celltype) & contains(Genotype)
                    ,CellTypes[1,],tad,"from.baitEnsemblID","to.oePromoterEnsemblID") %>%
      dplyr::filter(!is.na(to.oePromoterSymbol))%>%
      dplyr::rename(score=5,bait=3,oeProm=4)%>%
      dplyr::filter(bait !="negative" & oeProm!="negative")%>%
      dplyr::filter(!(from.baitSymbol==to.oePromoterSymbol))
    
    if(ALL==FALSE){Sample<-Sample%>%filter(score>=5)}
    
    links <- Sample %>%
      dplyr::select("from.baitSymbol","to.oePromoterSymbol",score,tad)%>%
      group_by(from.baitSymbol,to.oePromoterSymbol)%>%
      mutate(Score=mean(score))%>%
      mutate(score=Score)%>%
      distinct(.,from.baitSymbol,to.oePromoterSymbol,.keep_all=TRUE)
    
    Baits <- Sample %>% 
      dplyr::select("from.baitSymbol",bait,"from.baitEnsemblID") %>%
      dplyr::rename(Categories=2,Genes=1,EnsembleID=3)
    
    OeProm <- Sample %>%
      dplyr::select("to.oePromoterSymbol",oeProm,"to.oePromoterEnsemblID") %>%
      dplyr::rename(Categories=2,Genes=1,EnsembleID=3)
    
    Genes <- rbind(Baits,OeProm[!(OeProm$Genes%in%Baits$Genes),]) %>%
      dplyr::filter(!is.na(Categories))%>%
      dplyr::distinct(Genes,.keep_all=TRUE)
    
    if(AceGradient){Genes <- Genes%>%
      dplyr::left_join(Prom_df%>%tibble::rownames_to_column("EnsembleID"),by="EnsembleID")%>%
      dplyr::distinct(Genes,.keep_all=TRUE)
    }
    
    
    
    
  }else{ #PE Networks
    
    Sample <- cbind(rdat,score)%>%
      as.data.frame() %>%
      filter(type.oe=="enhancer")%>%
      dplyr::select("from.baitSymbol","category.oeEnh.AtacID",
                    starts_with("category.bait.") & contains(Time) & contains(Celltype) & contains(Genotype),
                    starts_with("category.oeEnh.") & contains(Time) & contains(Celltype) & contains(Genotype)
                    ,CellTypes[1,],tad,"from.baitEnsemblID","category.oeEnh.AtacCord") %>%
      dplyr::rename(score=5,Bait=3,Enh=4) %>%
      dplyr::filter(!is.na(Enh) & Enh != "negative" & !is.na(Bait))
    
    if(ALL==FALSE){Sample<-Sample%>%filter(score>=5)}
    
    links <- Sample %>%
      dplyr::select("from.baitSymbol","category.oeEnh.AtacID",score,tad)%>%
      group_by(from.baitSymbol,category.oeEnh.AtacID)%>%
      mutate(Score=mean(score))%>%
      mutate(score=Score)%>%
      distinct(.,from.baitSymbol,category.oeEnh.AtacID,.keep_all=TRUE)
    
    
    Baits <- Sample %>% 
      dplyr::select("from.baitSymbol",Bait,"from.baitEnsemblID") %>%
      dplyr::rename(Categories=2,Genes=1,EnsembleID=3) %>%
      dplyr::mutate(Type="P")
    
    
    if(AceGradient){Baits<-Baits%>%
      dplyr::left_join(Prom_df%>%tibble::rownames_to_column("EnsembleID"),by="EnsembleID")}     
    
    OeEnh <- Sample %>%
      dplyr::select("category.oeEnh.AtacID",Enh,"category.oeEnh.AtacCord") %>%
      dplyr::rename(Categories=2,Genes=1,EnsembleID=3)%>%
      dplyr::mutate(Type="E")
    
    
    if(AceGradient){OeEnh<-OeEnh%>%
      dplyr::left_join(Enh_df%>%tibble::rownames_to_column("EnsembleID"),by="EnsembleID")
    }
    
    Genes <- rbind(Baits,OeEnh) %>%
      dplyr::distinct(Genes,.keep_all=TRUE)
  }         
  Links =list(links,Genes)
  names(Links)=c("links","Genes")
  return(Links)
}

########
#Network Plot
#############
rescale = function(x,a,b,c,d){c + (x-a)/(b-a)*(d-c)}

NetworkPlot = function (data=SE,single=FALSE,gene="Hoxa2",Celltype="FNP",Genotype="WT",Time="E11.5",P_P=TRUE,size=FALSE,TAD=FALSE,ALL=FALSE,AceGradient=FALSE,Enh_size=4){
  if(Genotype=="WT"){CellTypes=subset(colData(data),celltype==Celltype & genotype==tolower(Genotype) & time==Time, select="sample")
  }else{CellTypes=subset(colData(data),celltype==Celltype & genotype==Genotype & time==Time, select="sample")}
  
  # Setting Nodes colors
  if(P_P){
    VexColorLabel=c("#fcac62","#66C2A5","#FC8D62","#8DA0CB","#ffd92f")
    names(VexColorLabel)=c("Pc-only","active","bivalent","negative","primed")  
  }else{
    VexColorLabel=c("#fcac62","#66C2A5","#FC8D62","#8DA0CB","#ffd92f","#FC8D62")
    names(VexColorLabel)=c("Pc-only","active","bivalent","negative","primed","poised")  
  }
  
  # Setting the range of size nodes
  if(dim(sub)[1]<500){Sizeranges=c(5,18)}else{Sizeranges=c(4,8)}
  
  liste =.NetworkTable(data,single,gene,Celltype,Genotype,Time,P_P,ALL=ALL,AceGradient=AceGradient)
  links=liste[["links"]]
  Genes=liste[["Genes"]]
  
  # if(P_P){ #remove self-loops
  #	links=as.data.frame(links[links[,1]!=links[,2],])
  # 	Genes=Genes[Genes$Genes%in%c(links[,1],links[,2]),]
  
  # }
  
  Genes$colors=VexColorLabel[Genes[,"Categories"]]
  ActiveGradient=alpha("#66C2A5",seq(0.2,1,length.out=5))
  names(ActiveGradient)=seq(1,5,length.out=5)
  if(AceGradient){
    for(i in 1:length(Genes$colors)){
      Genes$colors[i]=ifelse(Genes$Categories[i]=="active",ActiveGradient[Genes$quantile_H3K27ac[i]],Genes$colors[i])
    }
  }
  
  network <- graph_from_data_frame(d=links,directed=F)  
  Vcolors=Genes[V(network)$name%in%Genes$Genes,"colors"] 
  
  if(TAD){
    edge.color=ifelse(links$score>=5,ifelse(links$tad=="Intra","gray","purple"),"white")
    E(network)$color=edge.color
  }
  if(ALL){delete.edges(network, which(E(network)$score<5))}
  
  #Edge=ntile(links$score,6)
  Edge=scales::rescale(links$score,to=c(1,8))
  
  if(P_P==FALSE) {
    LabelColor=ifelse(Genes[,"Type"]=="P","darkblue",Vcolors)
  }else{LabelColor="darkblue"}
  
  if(single==FALSE){
    if(P_P & size) { 
      d=degree(network)
      #Vsize=rescale(d, min(d), max(d), 5, 15)
      Vsize=scales::rescale(d,to=Sizeranges)
      V(network)$label.cex=.4
    }else if(P_P & size==FALSE){
      Vsize=4
      V(network)$label.cex=.4
    }else if(P_P==FALSE & size==FALSE){
      Vsize=ifelse(Genes[,"Type"]=="P",Enh_size,Enh_size/2)
      V(network)$label.cex=ifelse(Genes[,"Type"]=="P",.4,.001)
    }
    
  }else{
    V(network)$label.cex=.6
    Vsize =40
    
  }
  
  
  Layout=layout.fruchterman.reingold(network)
  plot(network, vertex.size= Vsize ,vertex.label.color=LabelColor , vertex.color=Vcolors,vertex.frame.color = Vcolors, edge.width=Edge,edge.curved=0.5,vertex.label.family="Times",vertex.label.font=1,main=CellTypes[1,],layout=Layout)
  legend("bottomright" ,legend=names(VexColorLabel), col =VexColorLabel  , bty = "n", pch=20 , pt.cex = 1, cex = 0.8,text.col="black" , horiz = F)
  print(Genes[Genes$Genes%in%TFs2_symbol,]$Genes)
  if(P_P){
    Promoters = Legend(labels=names(VexColorLabel), legend_gp = gpar(fill=VexColorLabel), title = "Promoters Classes")
    AcethylationLevel= Legend(labels=c("Low","","Medium","","High"), legend_gp = gpar(fill=ActiveGradient), title = "Quantile H3K27ac") 
    pdf(paste0("19_legend_PP_",CellTypes[1,],".pdf"))
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nr = 1, nc = 2)))
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
    grid.draw(Promoters)
    upViewport()
    
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
    grid.draw(AcethylationLevel)
    upViewport()
    
    dev.off()
  }
  
  #print(links[links$from.baitSymbol%in%TF_symbol])
  return(Layout)
}


########
#Network Differential
#############


NetworkDiff= function(data=E11.5,gene="Hoxa2",Celltype="FNP",Genotype=c("WT","Ezh2KO"),Time="E11.5",P_P=FALSE) {
  single=TRUE
  VexColorLabel=c("#FF7F00","#66C2A5","#FC8D62","#8DA0CB","#FFD92F","#FC8D62")
  names(VexColorLabel)=c("Pc-only","active","bivalent","negative","primed","poised")  
  Vsize =40
  CellTypes=subset(colData(data),celltype%in%Celltype & genotype%in%Genotype & time%in%Time, select="sample")
  
  Diff=list()
  for ( i in 1:length(Genotype)) {
    Diff[[Genotype[i]]] =.NetworkTable(data,single,gene,Celltype,Genotype[i],Time,P_P)            
  }  
  
  #Union of Networks
  UniNet= graph_from_data_frame(d=rbind(Diff[["wt"]][["links"]],Diff[["Ezh2KO"]][["links"]]),directed=F)  
  
  #Get the coords of the nodes
  Coords <- layout_with_fr(UniNet) %>% 
    as_tibble %>%
    bind_cols(data_frame(names = names(V(UniNet))))
  
  #..Normalize the coords of nodes
  #coords_rescaled <- sapply(Coords[-3], function(x) -1+((x-min(x))*2)/diff(range(x))) 
  #rownames(coords_rescaled) <- Coords$names
  
  #..Construct each graph 
  LinksWt=Diff[["wt"]][["links"]]
  LinksEzh2ko=Diff[["Ezh2KO"]][["links"]]
  wtNet= graph_from_data_frame(d=LinksWt,directed=F)  
  ezh2koNet= graph_from_data_frame(d=LinksEzh2ko,directed=F)
  
  #..Get Color nodes for each Net
  GenesWt=Diff[["wt"]][["Genes"]]
  ColorsWt=VexColorLabel[GenesWt[,"Categories"]]
  
  GenesEzh2ko=Diff[["Ezh2KO"]][["Genes"]]
  ColorsEzh2ko=VexColorLabel[GenesEzh2ko[,"Categories"]]
  
  
  #..Set diff nodes +/- in ezh2ko compare to wt
  Plus = setdiff(GenesEzh2ko$Genes,GenesWt$Genes)
  Minus = setdiff(GenesWt$Genes,GenesEzh2ko$Genes)
  Intersect= intersect(GenesWt$Genes,GenesEzh2ko$Genes)
  
  #Reconstruct the Ezh2ko Net
  GenesEzh2ko$Frame = "black"
  LinksEzh2ko$Edge="grey"
  GenesEzh2ko[GenesEzh2ko$Genes %in% Plus,"Frame"]="red"
  LinksEzh2ko[LinksEzh2ko$category.oeEnh.AtacID %in% Plus,"Edge"]="red"
  
  #Reconstruct the Wt Net
  GenesWt$Frame = "black"
  LinksWt$Edge="grey"
  GenesWt[GenesWt$Genes %in% Minus,"Frame"]="darkgreen"
  LinksWt[LinksWt$category.oeEnh.AtacID %in% Minus,"Edge"]="darkgreen"
  
  
  par(mfrow=c(1,2),mar=c(2,2,2,2)) 
  
  wtNet %>%
    plot(., vertex.size= Vsize , vertex.color=ColorsWt,vertex.frame.color = GenesWt$Frame, edge.width=log2(LinksWt$score),edge.curved=0.5,edge.color=LinksWt$Edge,,vertex.label.family="Times",vertex.label.font=1,main=CellTypes[1,],layout=layout.fruchterman.reingold)
  
  ezh2koNet %>%
    plot(., vertex.size= Vsize , vertex.color=ColorsEzh2ko,vertex.frame.color = GenesEzh2ko$Frame, edge.width=log2(LinksEzh2ko$score),edge.curved=0.5,edge.color=LinksEzh2ko$Edge,vertex.label.family="Times",vertex.label.font=1,main=CellTypes[2,],layout=layout.fruchterman.reingold) 
}





##########
##Barplot 
#################

BarPlot= function(data=SE, gene="Hoxa2",celltypes="FNP",times="E10.5",genotypes="WT") {
  Total=data.frame()
  SubM <- data[rowData(data)$from.baitSymbol==gene,] 
  
  if(genotypes=="Ezh2cKO") {
    BarColLabel=c("#66C2A5","#8DA0CB","#FFD92F")
    names(BarColLabel)=c("active","negative","primed")  
    
    CellTypes=subset(colData(SubM),celltype%in%celltypes & genotype%in%genotypes & time%in%times, select="sample")
    for ( i in 1:dim(CellTypes)[1]) {
      Stats <- cbind(rowData(SubM) %>%
                       as.data.frame() %>%
                       dplyr::select(starts_with("category.oeEnh") & contains(!!CellTypes[i,])),assay(SubM)[,CellTypes[i,]]) %>%
        as.data.frame() %>%
        dplyr::rename(score=dim(.)[2],category.oeEnh=1) %>%
        dplyr::filter(score >=3 & !is.na(category.oeEnh)) %>%
        dplyr::group_by(category.oeEnh) %>%
        dplyr::summarise(Per=(n()*100)/dim(.)[1])%>%
        dplyr::mutate(cellType=CellTypes[i,])
      
      Total=rbind(Total,Stats)
    }
    ggplot(Total %>%
             dplyr::mutate(cellType= factor(cellType,levels=CellTypes[,"sample"])), aes(fill=category.oeEnh, y=Per, x=cellType)) + 
      geom_bar(position="stack", stat="identity") +  
      geom_col(fill=BarColLabel[Total$category.oeEnh]) + 
      scale_fill_manual(values =BarColLabel)+
      theme_cowplot() +
      coord_flip()
    
  }else{
    CellTypes=subset(colData(SubM),celltype%in%celltypes & genotype%in%tolower(genotypes) & time%in%times, select="sample")
    BarColLabel=c("#FF7F00","#66C2A5","#FC8D62","#8DA0CB","#FFD92F","#FC8D62")
    names(BarColLabel)=c("Pc-only","active","bivalent","negative","primed","poised")  
    
    for ( i in 1:dim(CellTypes)[1]) {
      Stats <- cbind(rowData(SubM) %>%
                       as.data.frame() %>%
                       dplyr::select(starts_with("category.oeEnh") & contains(!!CellTypes[i,])),assay(SubM)[,CellTypes[i,]]) %>%
        as.data.frame() %>%
        dplyr::rename(score=dim(.)[2],category.oeEnh=1) %>%
        dplyr::filter(score >=3 & !is.na(category.oeEnh)) %>%
        dplyr::group_by(category.oeEnh) %>%
        dplyr::summarise(Per=(n()*100)/dim(.)[1]) %>%
        dplyr::mutate(cellType=CellTypes[i,])  
      
      Total=rbind(Total,Stats)
    }
    
    ggplot(Total %>%
             dplyr::mutate(cellType= factor(cellType,levels=CellTypes[,"sample"])), aes(fill=category.oeEnh, y=Per, x=cellType)) + 
      geom_bar(position="stack", stat="identity") +  
      geom_col(fill=BarColLabel[Total$category.oeEnh]) + 
      scale_fill_manual(values =BarColLabel)+
      ggtitle(paste("% of different Enhancers of",gene,times,genotypes))+
      theme_cowplot() +
      coord_flip()
  }
}       


########
#Network similarity
#############

Similarity = function(data=E11.5,genes="Hoxa2",times="E11.5",genotypes="WT",celltypes="FNP") { 
  data <- data[rowData(data)$from.baitSymbol==genes,] 
  rdat <- rowData(data)
  score <- assay (data)
  Experiment <- colData(data) 
  CellTypes=subset(Experiment,celltype%in%celltypes & genotype%in%genotypes & time%in%times, select="sample")
  networkList=list()
  for ( i in 1:length(celltypes)) {
    if (genotypes=="Ezh2KO") {
      Sample <- cbind(rdat,score)%>%
        as.data.frame() %>%
        dplyr::select("from.baitSymbol","category.oeEnh.AtacID",
                      starts_with("category.bait.") & contains(celltypes[i]),
                      starts_with("category.oeEnh.") & contains(celltypes[i])
                      ,starts_with(celltypes[i]) & contains("KO")) %>%
        dplyr::rename(score=5,Bait=3,Enh=4) %>%
        dplyr::filter(score >= 3) %>%
        dplyr::filter(!is.na(Enh) & Enh != "negative")
      
    } else {
      Sample <- cbind(rdat,score)%>%
        as.data.frame() %>%
        dplyr::select("from.baitSymbol","category.oeEnh.AtacID",
                      starts_with("category.bait.") & contains(celltypes[i]),
                      starts_with("category.oeEnh.") & contains(celltypes[i])
                      ,starts_with(celltypes[i])) %>%
        dplyr::rename(score=5,Bait=3,Enh=4) %>%
        dplyr::filter(score >= 3) %>%
        dplyr::filter(!is.na(Enh) & Enh != "negative")
    }
    
    links <- Sample %>%
      dplyr::select("from.baitSymbol","category.oeEnh.AtacID",score)
    
    networkList[[celltypes[i]]] <- graph_from_data_frame(d=links,directed=F)  
  }
  
  E8.5 <- E10.5[rowData(E10.5)$from.baitSymbol==genes & assay(E10.5)[,"E8.5_progenitors"] >=3,]
  Proge <-  cbind(rowData(E8.5),assay(E8.5))%>%
    as.data.frame() %>%
    dplyr::select("from.baitSymbol","category.oeEnh.AtacID",
                  "category.bait.E8.5_progenitors","category.oeEnh.E8.5_progenitors"
                  ,"E8.5_progenitors") %>%
    dplyr::filter(!is.na(category.oeEnh.E8.5_progenitors) & category.oeEnh.E8.5_progenitors != "negative")
  
  links <- Proge %>%
    dplyr::select("from.baitSymbol","category.oeEnh.AtacID","E8.5_progenitors")
  
  networkList[["E8.5_progenitors"]] <- graph_from_data_frame(d=links,directed=F)           
  
  fracEdgesRowAlsoInCol <- communitySimilarity <- matrix(NA, nrow = length(Samples), ncol = 1,
                                                         dimnames = list(Samples,"E8.5_progenitors"))
  
  for (smpl1 in celltypes) {
    
    m1 <- components(networkList[[smpl1]])$membership
    
    m2 <- components(networkList[["E8.5_progenitors"]])$membership
    
    # compare graphs (fraction of edges in g1 also contain in g2)
    dg <- difference(networkList[[smpl1]], networkList[["E8.5_progenitors"]])
    fracEdgesRowAlsoInCol[smpl1, "E8.5_progenitors"] <- 1 - (length(E(dg)) / length(E(networkList[[smpl1]])))
    
    
    # ... find the common nodes
    cids <- intersect(names(m1), names(m2))
    communitySimilarity[smpl1, "E8.5_progenitors"] <- (length(cids)*100)/ length(names(m2))
    
  }
  
  ggplot(communitySimilarity %>%
           as.data.frame() %>%
           dplyr::mutate(CellType =factor(rownames(communitySimilarity),levels=c("FNP","Mx","Md","PA2")))
         , aes(y=E8.5_progenitors, x=CellType)) + 
    geom_bar(position="stack", stat="identity") +  
    xlab("Cell Types") +
    ylab("%") +
    ggtitle(paste("% of common nodes",genes,times,genotypes,"relative to progenitors"))+
    theme_cowplot() +
    coord_flip() 
}


########
#Expression Plot
#################

ExprPlot =function(data=E10.5,gene ="Hoxa2",time="E10.5",genotype="WT") {
  
  Expression=rowData(data) %>%
    as.data.frame() %>%
    dplyr::filter(from.baitSymbol==gene) %>%
    dplyr::select(starts_with("expression")) 
  
  if (time =="E10.5") {
    Expression <- Expression %>%
      dplyr::relocate(contains("E8.5")) 
  }
  if(time=="E11.5" & tolower(genotype)=="wt") {
    Expression <- Expression %>%
      dplyr::select(!contains("Ezh2cKO")) 
    
    
  }
  
  if(time=="E11.5" & genotype=="Ezh2cKO") {
    
    Expression <- Expression %>%
      dplyr::select(contains("Ezh2cKO")) 
    
  }
  
  colnames(Expression)= sub("expression.","",colnames(Expression))
  ColNames <-colnames(Expression) 
  Expression <- Expression%>%
    dplyr::distinct()%>%
    tidyr::pivot_longer(cols=ends_with("WT"),names_to="CellType",values_to="Rna")%>%
    mutate(CellType=factor(CellType, levels=CellType))
  
  P<- ggplot(Expression, aes(x=CellType, y=Rna)) +
    geom_segment( aes(xend=CellType, yend=0),linewidth=10) +
    geom_point( size=12, color="purple") +
    coord_flip() +
    theme_bw() +
    ggtitle(gene)+
    ylab("Mean(Rpkm)")+
    xlab("Cell-Type")+
    theme(axis.text=element_text(size=20,face="bold"),
          axis.title=element_text(size=16),
          plot.title = element_text(color="red", size=20, face="bold.italic"))
  
  
  return(P)
}
