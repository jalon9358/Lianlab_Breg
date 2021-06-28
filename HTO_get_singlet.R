###############################################################
############      citeseq_divide_HTOtag       #################
###############################################################
library(Seurat)
library(magrittr)

##function part----
HTO_get_singlet=function(path,filename="",positive.quantile=0.99,width=6.5,height=5,kfunc="clara"){
  #import cellranger exported data and create seurat object
  data=Read10X(data.dir = path)  
  object=CreateSeuratObject(counts = data$`Gene Expression`)
  
  #cell-hash classification
  object[["Antibody"]]=CreateAssayObject(counts = data$`Antibody Capture`)
  object=NormalizeData(object,assay="Antibody",normalization.method = "CLR")
  object=HTODemux(object,assay = "Antibody",positive.quantile = positive.quantile,kfunc = kfunc)
  fre=as.data.frame(table(object$hash.ID)/ncol(object))
  count=as.data.frame(table(object$hash.ID))
  merge=merge(fre,count,all=T,by="Var1")
  colnames(merge)=c("Class","Frequency","Cell_count")
  write.csv(merge,paste0(filename,"_Antibody_classification.csv"),row.names = F)
  
  Idents(object)="hash.ID"
  pdf(paste0(filename,"_HTO_Ridgeplot.pdf"), width = width, height = height)
  p=RidgePlot(object, assay = "Antibody", features = rownames(object[["Antibody"]]), ncol = 1)
  print(p)
  dev.off()
  pdf(paste0(filename,"_HTO_scatter.pdf"), width = width, height = height)
  p=suppressWarnings(FeatureScatter(object, feature1 = "Sampletag1", feature2 = "Sampletag2"))
  print(p)
  dev.off()
  pdf(paste0(filename,"_HTO_violin.pdf"), width = width, height = height)
  p=VlnPlot(object, features = "nFeature_RNA", pt.size = 0.1, log = TRUE)
  print(p)
  dev.off()
  
  #remove doublet and negative cells 
  Idents(object)="Antibody_classification.global"
  object=subset(object, idents = "Singlet")
  return(object)
}

