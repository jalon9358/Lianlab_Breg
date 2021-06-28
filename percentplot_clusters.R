percentplot_clusters=function(seurat_object,filename,width=6,height=4,flip=F){
  num <- table(seurat_object$b10)
  cluster_name=levels(seurat_object)
  #cluster_name=sort(cluster_name)
  frequency_matrix=data.frame(names(num))
  group=frequency_matrix$names.num.
  for (i in cluster_name) {
    cell_name=row.names(seurat_object@meta.data[seurat_object@active.ident==i,])
    subset_data=subset(seurat_object, cells = cell_name)
    frequency=table(subset_data@meta.data$b10)/ncol(subset_data)
    frequency=as.data.frame(frequency)
    zero=setdiff(as.character(unique(seurat_object$b10)),
                 names(table(subset_data$b10)))
    zero=as.data.frame(zero)
    zero$Freq=c(rep(0,nrow(zero)))
    names(zero)[1]="Var1"
    frequency=rbind(frequency,zero)
    suppressWarnings(if(!is.na(as.numeric(as.character(frequency$Var1)[1]))){
      print("reorder cluster name")
      frequency=frequency[order(as.numeric(as.character(frequency$Var1))),]
    })
    row.names(frequency)=frequency$Var1
    frequency=frequency[group,]
    frequency_matrix[,i]=frequency$Freq
  }
  cluster_name2=as.character(cluster_name)
  names(frequency_matrix)=c("cluster_name",cluster_name2)
  row.names(frequency_matrix)=frequency_matrix[,1]
  write.csv(frequency_matrix,paste0(filename,"_frequency_matrix.csv"),row.names = F)
  library(reshape2)
  data_rownames <- row.names(frequency_matrix)
  data_colnames <- colnames(frequency_matrix)
  frequency_matrix$cluster <- data_rownames
  frequency_matrix=frequency_matrix[,-1]
  frequency_matrix$cluster=factor(frequency_matrix$cluster,levels = frequency_matrix$cluster)
  frequency_matrix_m <- reshape2::melt(frequency_matrix, id.vars=c("cluster"))
  frequency_matrix_m <- frequency_matrix_m %>% dplyr::group_by(variable) %>% 
    dplyr::mutate(count=sum(value)) %>% dplyr::mutate(freq=round(100*value/count,2))
  #frequency_matrix_m=frequency_matrix_m[sort(frequency_matrix_m$Clusters),]
  colnames(frequency_matrix_m)=c("Group","Cluster","Frequency","count","Percent")
  colorm=cbind(names(num),color_blue[5:length(names(num))])
  row.names(colorm)=colorm[,1]
  colorm=colorm[sort(names(num)),]
  color_percent=as.character(colorm[,2])
  color_percent=colors2[2:1]
  p <- ggplot(frequency_matrix_m, aes(x=Cluster, y=Percent, group=Group)) +
    geom_bar(stat="identity", color="white",position="fill", aes(fill=Group))+
    scale_fill_manual(values=color_percent,breaks=data_rownames)+theme_bw()+
    theme(panel.border = element_blank(),panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
          axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))
  if(flip==TRUE){
    p=p+coord_flip()
  }
  pdf(paste0("percent_plot_",filename,".pdf"),width =width, height = height)
  print(p)
  dev.off()
}