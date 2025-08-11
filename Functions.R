suppressPackageStartupMessages({
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

NetworkPlot = function (data=SE,single=FALSE,gene="Hoxa2",Celltype="FNP",Genotype="WT",Time="E11.5",P_P=TRUE,size=FALSE,Enh_size=4){
 
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
    
  Edge=scales::rescale(links$score,to=c(1,8))
  
  if(P_P==FALSE) {
    LabelColor=ifelse(Genes[,"Type"]=="P","darkblue",Vcolors)
  }else{LabelColor="darkblue"}
  
  if(single==FALSE){
    if(P_P & size) { 
      d=degree(network)
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

}



