
clonotype_share_samples=function(seurat_object){
  BCR_matrix=subset(seurat_object@meta.data,subset=(clonotype!="None"))
  
  BCR_matrix_split=split(BCR_matrix,BCR_matrix$samples)
  cluster_no=length(unique(BCR_matrix$samples))
  count_shared=as.data.frame(matrix(c(rep(0,cluster_no^2)),ncol = cluster_no,nrow = cluster_no))
  rownames(count_shared)=names(BCR_matrix_split)
  colnames(count_shared)=names(BCR_matrix_split)
  clonotype_shared=count_shared
  
  
  for (i in 1:length(BCR_matrix_split)) {
    for (j in 1:length(BCR_matrix_split)){
      inter=intersect(BCR_matrix_split[[i]]$clonotype,BCR_matrix_split[[j]]$clonotype)
      clonotype_shared[i,j]=length(inter)
      count=0
      for (k in inter) {
        count=count+min(length(which(BCR_matrix_split[[i]]$clonotype==k)),
                        length(which(BCR_matrix_split[[j]]$clonotype==k)))
      }
      count_shared[i,j]=count
    }
  }
  
  for (i in 1:length(BCR_matrix_split)) {
    count_shared[i,i]=0
    clonotype_shared[i,i]=0
  }
  
  p1=pheatmap(count_shared,display_numbers = T,number_format = "%.0f",number_color = "black",
              fontsize_number = 20,cluster_rows = F,cluster_cols = F,border_color = "black",
              color = colorRampPalette(colors = c("white","#0076C5","#F99400","firebrick"))(100))
  pdf("count_shared_samples.pdf",width=6,height = 5)
  print(p1)
  dev.off()
  p2=pheatmap(clonotype_shared,display_numbers = T,number_format = "%.0f",number_color = "black",
              fontsize_number = 20,cluster_rows = F,cluster_cols = F,border_color = "black",
              color = colorRampPalette(colors = c("white","#0076C5","#F99400","firebrick"))(100))
  pdf("clonotype_shared_samples.pdf",width=6,height = 5)
  print(p2)
  dev.off()
}

clonotype_share_clusters=function(seurat_object){
  BCR_matrix=subset(seurat_object@meta.data,subset=(clonotype!="None"))
  
  BCR_matrix_split=split(BCR_matrix,BCR_matrix$seurat_clusters)
  cluster_no=length(unique(BCR_matrix$seurat_clusters))
  count_shared=as.data.frame(matrix(c(rep(0,cluster_no^2)),ncol = cluster_no,nrow = cluster_no))
  rownames(count_shared)=names(BCR_matrix_split)
  colnames(count_shared)=names(BCR_matrix_split)
  clonotype_shared=count_shared
  
  for (i in 1:length(BCR_matrix_split)) {
    for (j in 1:length(BCR_matrix_split)){
      inter=intersect(BCR_matrix_split[[i]]$clonotype,BCR_matrix_split[[j]]$clonotype)
      clonotype_shared[i,j]=length(inter)
      count=0
      for (k in inter) {
        count=count+min(length(which(BCR_matrix_split[[i]]$clonotype==k)),
                        length(which(BCR_matrix_split[[j]]$clonotype==k)))
      }
      count_shared[i,j]=count
    }
  }
  
  for (i in 1:length(BCR_matrix_split)) {
    count_shared[i,i]=0
    clonotype_shared[i,i]=0
  }
  
  p1=pheatmap(count_shared,display_numbers = T,number_format = "%.0f",number_color = "black",
              fontsize_number = 20,cluster_rows = F,cluster_cols = F,border_color = "black",
              color = colorRampPalette(colors = c("white","#0076C5","#F99400","firebrick"))(100))
  pdf("count_shared_clusters.pdf",width=6,height = 5)
  print(p1)
  dev.off()
  p2=pheatmap(clonotype_shared,display_numbers = T,number_format = "%.0f",number_color = "black",
              fontsize_number = 20,cluster_rows = F,cluster_cols = F,border_color = "black",
              color = colorRampPalette(colors = c("white","#0076C5","#F99400","firebrick"))(100))
  pdf("clonotype_shared_clusters.pdf",width=6,height = 5)
  print(p2)
  dev.off()
}


