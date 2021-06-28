###source funcions and packages----
source("/Users/jalon/Desktop/newest scripts/functions/source.R")

###import 10x scRNAseq data with HTO hashtag,here we use HTODemux funtion to divide different organs##
###and remove the doublets which labeled with two different hashtags###
path1="/Users/jalon/Desktop/help_analysis/YSY_B/20201209/S1/outs/count/filtered_feature_bc_matrix"
path2="/Users/jalon/Desktop/help_analysis/YSY_B/20201209/S2/outs/count/filtered_feature_bc_matrix"
path3="/Users/jalon/Desktop/help_analysis/YSY_B/20201209/S3/outs/count/filtered_feature_bc_matrix"
s1.singlet=HTO_get_singlet(path = path1,filename = "S1",positive.quantile = 0.95)
s2.singlet=HTO_get_singlet(path = path2,filename = "S2",positive.quantile = 0.96)
s3.singlet=HTO_get_singlet(path = path3,filename = "S3",positive.quantile = 0.96)

###add samples information
s1.singlet$samples=paste0("S1_",s1.singlet$Antibody_classification)
s2.singlet$samples=paste0("S2_",s2.singlet$Antibody_classification)
s3.singlet$samples=paste0("S3_",s3.singlet$Antibody_classification)

#merge all objects and add organs information
all.singlet <- merge(s1.singlet, y = c(s2.singlet,s3.singlet), 
                     add.cell.ids = c("S1", "S2","S3"), project = "Bcells")
all.singlet$samples=ifelse(all.singlet$samples=="S1_Sampletag1","Liver",
                           ifelse(all.singlet$samples=="S1_Sampletag2","Spleen",
                                  ifelse(all.singlet$samples=="S2_Sampletag1","mLN",
                                                ifelse(all.singlet$samples=="S3_Sampletag1","BM","PC"))))
#calculate mitochodrial and ribosome related gene frequency
all.singlet[["percent.mt"]] <- PercentageFeatureSet(all.singlet, pattern = "^mt-")
all.singlet[["percent.ribo"]] <- PercentageFeatureSet(all.singlet, pattern = "^Rp[s|l]")
#calculate log10GenesPerUMI
all.singlet$log10GenesPerUMI <- log10(all.singlet$nFeature_RNA)/log10(all.singlet$nCount_RNA)

##set the width and height of plots
width.ppi=6.5
height.ppi=5

##Figure S1B-plot the log10GenesPerUMI----
pdf("log10GenesPerUMI_QC.pdf", width = width.ppi*1.1, height = height.ppi)
p=ggplot(all.singlet@meta.data,aes(x=log10GenesPerUMI, color = samples, fill=samples)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8)+ggtitle(label = "Percent 99.8%")
print(p)
dev.off()

##data filter,choose the cells with nFeature_RNA more than 300 but less than 5000,
##log10GenesPerUMI more than 0.8 and percent.mt more than 10
all.singlet <- subset(all.singlet, subset=nFeature_RNA>300 & nFeature_RNA<5000 &
                         log10GenesPerUMI>0.8 & percent.mt<10 ) 

#filter genes,only choose the genes which expressed in more than 10 cells
counts=GetAssayData(object = all.singlet, slot = "counts")
all.singlet=CreateSeuratObject(counts,min.cells = 10,meta.data = all.singlet@meta.data)

#import BCR data
BCR.aggr=read.csv("/Users/jalon/Desktop/help_analysis/YSY_B/20201209/aggr0109/vdj_b/filtered_contig_annotations.csv")
#change cell barcode names which match the cellname in all.singlet seurat object 
BCR.aggr$barcode=paste0(ifelse(BCR.aggr$origin=="Liver_spleen","S1_",
                               ifelse(BCR.aggr$origin=="mLN","S2_","S3_")),BCR.aggr$barcode)
BCR.aggr$barcode=paste0(str_split_fixed(BCR.aggr$barcode,pattern = "-",2)[,1],"-1")
#classify the cells according to the chains number 
chain4=subset(BCR.aggr,subset=barcode%in%names(which(table(BCR.aggr$barcode)==4)))
chain1=subset(BCR.aggr,subset=barcode%in%names(which(table(BCR.aggr$barcode)==1)))
chain3=subset(BCR.aggr,subset=barcode%in%names(which(table(BCR.aggr$barcode)==3)))
chain2=subset(BCR.aggr,subset=barcode%in%names(which(table(BCR.aggr$barcode)==2)))
chain4_cellname=names(table(chain4$barcode))
chain1_cellname=names(table(chain1$barcode))
chain3_cellname=names(table(chain3$barcode))
chain2_cellname=names(table(chain2$barcode))
#insert BCR chains number information into seurat object
all.singlet$BCR_nchain=ifelse(colnames(all.singlet)%in%chain1cellname,"c1",
                              ifelse(colnames(all.singlet)%in%chain3_cellname,"c3",
                                     ifelse(colnames(all.singlet)%in%chain4_cellname,"c4",
                                            ifelse(colnames(all.singlet)%in%chain2_cellname,"c2","c0"))))
all.singlet$BCR_nchain2=ifelse(all.singlet$BCR_nchain=="c4","c4","Others")
Idents(all.singlet)="BCR_nchain"
##Figure S1C-violin plot the nFeature_RNA and nCount_RNA of B cells with c4 or others BCR chains----
pdf("Vlnplot_QC_BCR_c4.pdf", width = width.ppi*1.8, height = height.ppi)
VlnPlot(all.singlet, features = c("nFeature_RNA", "nCount_RNA"), 
        ncol = 2,pt.size = 0,cols = colors2)
dev.off()

#deplete the doublets with 4 BCR chains which are composed of two different heavy chains and two different light chains
all.singlet=subset(all.singlet,subset=BCR_nchain!="c4")

#insert BCR clonotypes information into seurat object
BCR.aggr.clono=subset(BCR.aggr,select=c(barcode,raw_clonotype_id))
BCR.aggr.clono=BCR.aggr.clono[!duplicated(BCR.aggr.clono$barcode),]
row.names(BCR.aggr.clono)=BCR.aggr.clono$barcode
meta=all.singlet@meta.data
meta$barcode=row.names(meta)
merge_barcode=merge(meta,BCR.aggr.clono,by="barcode",all.x=TRUE)
merge_barcode$raw_clonotype_id[which(is.na(merge_barcode$raw_clonotype_id))]="None"
all.singlet$clonotype=merge_barcode$raw_clonotype_id

#set the samples information as active idents 
Idents(all.singlet)="samples"

#Figure S1D-plot nFeature_RNA and percent.mt of different organs----
pdf("Vlnplot_QC_after.pdf", width = width.ppi, height = height.ppi*2)
p=VlnPlot(all.singlet, features = c("nFeature_RNA", "percent.mt"), ncol = 1,pt.size = 0,cols = colors2)
print(p)
dev.off()

cellnumber=table(all_singlet$samples)
cellnumber_mtx=data.frame(organs=names(cellnumber),cellcount=as.numeric(cellnumber))
##Figure S1E-barplot of the cell counts of different organs----
pdf("cellcount_organs.pdf",width = width.ppi,height = height.ppi)
p=ggplot(cellnumber_mtx,aes(x=organs,y=cellcount,fill=organs))+geom_bar(stat = "identity")+
  scale_fill_manual(values=colors2)+theme_bw()+
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"),
        axis.text.x = element_text( vjust = 0.5, hjust = 0.5, angle = 45))+
  geom_text(aes(label=cellcount),vjust=-0.5)
print(p)
dev.off()

#dimension reduction and clustering
all.singlet=all.singlet%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%
  RunPCA()
ElbowPlot(object = all.singlet,ndims = 50)
all.singlet <- FindNeighbors(object = all.singlet, dims = 1:20)
all.singlet <- FindClusters(object = all.singlet, resolution = 0.4)
all.singlet <- RunTSNE(object = all.singlet, dims = 1:20)
all.singlet <- RunUMAP(all.singlet, reduction = "pca", dims = 1:20)

#import breg gene sets
Breg_gene=read.table("/Users/jalon/Desktop/help_analysis/YSY_B/20210223/Breg_gene.txt",
                          sep = "\t",header = T)
Breg_gene=Breg_gene$Symbol

#calculate the gene signature score of Breg
#Figure 1C-export breg cell score plot----
breg_score=mymodule_score(seurat_object = all.singlet,
                          genelist = Breg_gene,filename = "Breg",
                          cutoff = 0.16,width.ppi =5 )

#insert breg label into all merged data
all.singlet$breg=ifelse(colnames(all.singlet)%in%colnames(Breg.cut),"Breg","Non-Breg")
#set breg label as the active ident
Idents(all.singlet)="breg"

#Figure 1D-plot Breg cells in total B cells----
pdf("Breg in total.pdf",width = 6.5,height = 5)
DimPlot(all.singlet,reduction = "tsne",group.by = "cluster",
        cols = c("hotpink2","lightgrey"))
dev.off()
#Figure 1E-split Breg cells plot by organs----
pdf("Breg in total samples.pdf",width = width.ppi*2.6,height = height.ppi*2)
DimPlot(all.singlet,reduction = "tsne",group.by = "cluster",
        cols = c(colors2[c(rep(1,7))],"lightgrey"),split.by = "samples",ncol = 3,pt.size = 1)
dev.off()
#build new label by merge samples and breg labels
all.singlet$organ_breg=paste0(all.singlet$samples,"_",all.singlet$breg)
Idents(all.singlet)="organ_breg"

#calculate upregulated genes in breg cells vs non-breg cells in various organs
sb_gene=FindMarkers(all.singlet,ident.1 = "Spleen_Breg",ident.2 = "Spleen_Non-Breg",
                    logfc.threshold = 0.5)
sb_upgene=rownames(subset(sb_gene,subset=avg_log2FC>0&p_val_adj<0.05))

PC_gene=FindMarkers(all.singlet,ident.1 = "PC_Breg",ident.2 = "PC_Non-Breg",
                    logfc.threshold = 0.5)
PC_upgene=rownames(subset(PC_gene,subset=avg_log2FC>0&p_val_adj<0.05))

l_gene=FindMarkers(all.singlet,ident.1 = "Liver_Breg",ident.2 = "Liver_Non-Breg",
                   logfc.threshold = 0.5)
l_upgene=rownames(subset(l_gene,subset=avg_log2FC>0&p_val_adj<0.05))

mLN_gene=FindMarkers(all.singlet,ident.1 = "mLN_Breg",ident.2 = "mLN_Non-Breg",
                     logfc.threshold = 0.5)
mLN_upgene=rownames(subset(mLN_gene,subset=avg_log2FC>0&p_val_adj<0.05))

BM_gene=FindMarkers(all.singlet,ident.1 = "BM_Breg",ident.2 = "BM_Non-Breg",
                    logfc.threshold = 0.5)
BM_upgene=rownames(subset(BM_gene,subset=avg_log2FC>0&p_val_adj<0.05))

#calculate the overlapped up-regulated genes in Breg cells in different organs
inter_breg_gene =intersect(BM_upgene,intersect(l_upgene,intersect(mLN_upgene,
                                                             intersect(l_upgene,intersect(sb_upgene,PC_upgene)))))
library(VennDiagram)
#Figure 2A-venn diagram of up-regulated genes in breg cells of various organs----
venn.diagram(list(Liver=l_upgene,BM=BM_upgene,mLN=mLN_upgene,Spleen=sb_upgene,PC=PC_upgene), 
             fill=colors[1:5], 
             alpha=c(0.5,0.5,0.5,0.5,0.5), cex=1, cat.fontface=4, filename="breg_gene_venn.png",
             imagetype = "png")

###Figure 2B-breg markers violin plot----
levels(all.singlet)=c("Spleen_Non-Breg","Spleen_Breg","Liver_Non-Breg","Liver_Breg","mLN_Non-Breg",   
                     "mLN_Breg", "PC_Non-Breg" ,"PC_Breg", "BM_Non-Breg"  , "BM_Breg" )
pdf("Breg_marker.pdf",width = 12,height = 9)
VlnPlot(all.singlet,features = c("Fcrl5","Cd9","Zbtb20","Ccnd2","Ccdc28b","Ptpn22"),pt.size = 0,cols = colors2,
        ncol = 2)
dev.off()

##SCENIC infer transcription network
library(SCENIC)
Idents(all.singlet)="breg"
all.singlet$seurat_clusters=all.singlet@active.ident
all.singlet$seurat_clusters= factor(all.singlet$seurat_clusters,levels = levels(Idents(all.singlet))) 
#randomly select 500 cells in breg cluster and non-breg cluster
nonbreg_name=row.names(subset(all.singlet@meta.data,
                         subset=seurat_clusters=="Non-Breg")[sample(nrow(subset(all.singlet@meta.data,
                                                                                subset=seurat_clusters=="Non-Breg")),500),])
bregcell_name=row.names(subset(all.singlet@meta.data,
                         subset=seurat_clusters=="Breg")[sample(nrow(subset(all.singlet@meta.data,
                                                                            subset=seurat_clusters=="Breg")),500),])
tran_Breg_cellname=c(nonbreg_name=,bregcell_name)
tran_Breg=subset(all.singlet,cells=tran_Breg_cellname)
exprMat  <-  as.matrix(tran_Breg@assays$RNA@counts)
tran_Breg$clusters=tran_Breg$seurat_clusters
col_use=c("samples","seurat_clusters","clusters")
cellInfo <-  tran_Breg@meta.data[,col_use]
colnames(cellInfo)=c('sample', 'cluster' ,'celltype')
saveRDS(cellInfo, file="int/cellInfo.Rds")
scenicOptions <- initializeScenic(org="mgi", 
                                  dbDir="/Users/jalon/Desktop/help_analysis/YSY_B/analysis_1210/mm9_cisTarget", 
                                  nCores=12) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 
### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions)
exprMat_filtered <- exprMat[genesKept, ]
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered_log, scenicOptions)

exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings

scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget"))
library(doParallel)
library(AUCell)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC")
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
# Read information from loom file:
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$cluster),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))

##Figure 2C-heatmap of upregulated transcription regulatory internet in Breg cells---- 
pdf("transcription_net_scenic.pdf",width = width.ppi*1,height = height.ppi*2)
ComplexHeatmap::Heatmap(regulonActivity_byCellType_Scaled, name="Regulon activity")
dev.off()

#Figure 2D-feature plot of Atf3, Myc and Sox5 genes----
myfeatureplot(all.singlet,genename = c("Atf3","Myc","Sox5"),filename = "tran")

breg_tran=all.singlet
Idents(breg_tran)="breg"
breg_tran$seurat_clusters=breg_tran$breg
#select Atf3 related genes
tran_gene=intersect(regulons$Atf3,rownames(breg_tran))
breg_tran=ScaleData(breg_tran)
#Figure 2E-heatmap of Atf3 related genes in Breg cells and non-Breg cells----
pdf("tran_atf3_3cluster.pdf",width =8, height = 6)
f1=DoHeatmap(breg_tran, features = tran_gene,
             label = TRUE,size=3,assay = "RNA",lines.width = 1,
             group.colors = colors,draw.lines = T)+
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))
print(f1)
dev.off()

#subset Breg cells 
Breg.cut=subset(breg.object,subset=score_cluster=="Positive")
Breg.counts=GetAssayData(object = Breg.cut, slot = "counts")
Breg.cut=CreateSeuratObject(Breg.counts,min.cells = 3,meta.data = Breg.cut@meta.data)
#dimension reduction and clustering
Breg.cut=Breg.cut%>%NormalizeData()%>%ScaleData()%>%FindVariableFeatures()%>%
  RunPCA()
ElbowPlot(object = Breg.cut,ndims = 50)
Breg.cut <- FindNeighbors(object = Breg.cut, dims = 1:30)
Breg.cut <- FindClusters(object = Breg.cut, resolution = 0.5)
Breg.cut <- RunTSNE(object = Breg.cut, dims = 1:30)
Breg.cut <- RunUMAP(Breg.cut, reduction = "pca", dims = 1:30)
##Figure 3A-tsne plot of Breg cell clusters----
pdf("Bregcut_plot_tsne_nolable.pdf",width = width.ppi*1,height = height.ppi)
DimPlot(object = Breg.cut, reduction = 'tsne',label = F,cols = colors2,pt.size = 1)
dev.off()
##Figure 3B-heatmap of top50 genes in each breg cell clusters---- 
myfindmarkers(Breg.cut,filename = "Breg",species = "mouse")
##Figure 3C-violin plot of selected masker genes of breg cell clusters----
violin_plot_gene=c("Ccr7","Apoe","Ighv11-2","S100a4","Ighv1-55","Cr2","Stmn1")
stacked_violin_plot(gene = violin_plot_gene,seurat_object = Breg.cut,flip = F,Mean = F,col = colors2,
                    filename = "breg",width = 10)
##Figure 3D-heatmap of GSVA score of selected functional genesets----
#gene set pathway
gspath="/Users/jalon/Desktop/help_analysis/YSY_B/20210304/Breg_function.gmt"
pheatmap_score=GSVA_cluster_diff(seurat_object = Breg.cut,gspath = gspath,species = "mouse",
                                     supercell_no = 10,topn = 20,logFCcutoff = 0.01,adjPvalueCutoff = 0.2,filename = "new")

pdf(paste0("GSVA_cluster_diff_","breg",".pdf"),width = 1*8,height = 6)
p=pheatmap(pheatmap_score,fontsize_row = 6,cluster_rows = F,cluster_cols = F,
           treeheight_row = 0,
           color = rev(RColorBrewer::brewer.pal(n=11,name = "RdBu")))
print(p)
dev.off()

##Figure 3E-split the tsne plot of breg cell clusters by samples label(organs)----
pdf("Bregcut_plot_tsne_samples_nolable.pdf",width = width.ppi*2.6,height = height.ppi*2)
DimPlot(object = Breg.cut, reduction = 'tsne',label = F,cols = colors2,pt.size = 2,split.by = "samples",
        ncol = 3)
dev.off()
##Figure 3F-plot the percentage of clusters in every organs----
percentplot_sample(Breg.cut,filename = "breg",width = 8)

##BCR analysis----
#order the clonotype counts in Breg cells
sort(table(Breg.cut$clonotype),decreasing = T)
#labeled top10 BCR clonotypes in all Breg cells
Breg.cut$BCR_top10=ifelse(Breg.cut$clonotype=="clonotype1","C1",
                          ifelse(Breg.cut$clonotype=="clonotype2","C2",
                                 ifelse(Breg.cut$clonotype=="clonotype4","C3",
                                        ifelse(Breg.cut$clonotype=="clonotype3","C4",
                                               ifelse(Breg.cut$clonotype=="clonotype11","C5",
                                                      ifelse(Breg.cut$clonotype=="clonotype7","C6",
                                                             ifelse(Breg.cut$clonotype=="clonotype17","C7",
                                                                    ifelse(Breg.cut$clonotype=="clonotype13","C8",
                                                                           ifelse(Breg.cut$clonotype=="clonotype15","C9",
                                                                                  ifelse(Breg.cut$clonotype=="clonotype9","C10",
                                                                                         "Others"))))))))))
Breg.cut$BCR_top10=factor(Breg.cut$BCR_top10,levels = c(paste0("C",1:10),"Others"))
##Figure 4A-plot top10 BCR clonotypes in Breg cells tsneplot----
pdf("Breg_BCR_top10.pdf",width = 6.5,height = 5)
DimPlot(Breg.cut,reduction = "tsne",group.by = "BCR_top10",
        cols = c(colors2[1:10],"lightgrey"),pt.size = 1)
dev.off()
##Figure 4B-split top10 BCR clonotypes in Breg cells tsneplot by samples----
pdf("Breg_BCR_top10_samples.pdf",width = 6.5*2.5,height = 5*2)
DimPlot(Breg.cut,reduction = "tsne",group.by = "BCR_top10",
        cols = c(colors2[1:10],"lightgrey"),pt.size = 1,split.by = "samples")
dev.off()
##Figure 4C and 4D-heatmap of the shared clonotypes between different samples----
clonotype_share_samples(Breg.cut)

##STARTRAC analysis of BCR
library(Startrac)
library(tictoc)
Breg_bcr=Breg.cut@meta.data
Breg_bcr=subset(Breg_bcr,select = c(seurat_clusters,samples,clonotype))
Breg_bcr$patient=rep("mouse",nrow(Breg_bcr))
Breg_bcr$Cell_Name=row.names(Breg_bcr)
colnames(Breg_bcr)=c("majorCluster","loc","clone.id","patient","Cell_Name")
Breg_bcr=subset(Breg_bcr,subset = clone.id!="None")
tic("Startrac.run")
out=Startrac.run(Breg_bcr,proj = "Breg",verbose = F,n.perm = 10)
toc()

##Figure 4E-heatmap of the transition index score of different Breg cell clusters----
pdf("startrac_BCR_share.pdf",width = width.ppi,height = height.ppi)
plot(out,index.type="pairwise.tran",byPatient=F)
dev.off()
##Figure 4F and 4G-plot the expansion and migration index score of different Breg cell clusters----
pdf("Breg_startrac.pdf",width = width.ppi*0.9,height = height.ppi*1.8)
plot(out,index.type="cluster.all",byPatient=F)+scale_fill_manual(values = colors)
dev.off()

#import b10 gene set
B10_gene=read.table("/Users/jalon/Desktop/help_analysis/YSY_B/20210223/B10_gene.txt",
                     sep = "\t",header = T)
B10_gene=B10_gene$Symbol
#calculate B10 signature score in Breg cells
b10.score=mymodule_score(seurat_object = Breg.cut,genelist = b10_gene,cutoff = 0.35,width.ppi = 5,
                           filename = "b10")
breg.cellname=colnames(subset(breg.object,subset=score_cluster=="Positive"))
b10.cellname=colnames(subset(b10.score,subset=score_cluster=="Positive"))

##Figure 5A-Venn diagram of B10 cells and Breg cells----
venn.diagram(list(Breg=breg.cellname,B10=b10.cellname), fill=c("red","blue"), 
             alpha=c(0.5,0.5), cex=0, cat.fontface=4, filename="breg_venn_nonumber.png",
             imagetype = "png")
#add b10 label in Breg object
Breg.cut$b10=ifelse(colnames(Breg.cut)%in%b10.cellname,"B10","Non-B10")
##Figure 5B-plot B10 and non-B10 cells in all Breg cells tsneplot----
pdf("B10 in Breg.pdf",width = 6.5,height = 6)
DimPlot(Breg.cut,reduction = "tsne",group.by = "b10",
        cols = c(colors2[2:1]))
dev.off()

getwd()
##Figure 5C-barplot of the percentage of B10 and Non-B10 cells in different Breg clusters----
percentplot_clusters(seurat_object=Breg.cut,filename="b10_cluster_percent",width=6,height=4,flip=F)

##Figure 5D-pie plot of B10,non-B10 Breg and non-Breg cell percentages in different organs----
source("pie_b10.R")
pie_b10(seurat.object = all.singlet,filename = "b10_spleen_label",organ = "Spleen",label = T)
pie_b10(seurat.object = all.singlet,filename = "b10_BM_label",organ = "BM",label = T)
pie_b10(seurat.object = all.singlet,filename = "b10_PC_label",organ = "PC",label = T)
pie_b10(seurat.object = all.singlet,filename = "b10_Liver_label",organ = "Liver",label = T)
pie_b10(seurat.object = all.singlet,filename = "b10_mLN_label",organ = "mLN",label = T)

##Figure 5E-volcano show differential expression genes between B10 and non-B10 Breg cells----
Breg_only_marker=FindMarkers(Breg.cut,ident.1 = "Non-B10")
myvolcano(as.data.frame(Breg_only_marker),logfc.cutoff=0.5,p.cutoff=0.01,filename="breg_b10",
          width=6.5,height=5,gene.plot = 10,pt.size = 2)
myvolcano(as.data.frame(Breg_only_marker),logfc.cutoff=0.5,p.cutoff=0.01,filename="breg_b10_nolabel",
          width=6.5,height=5,gene.plot = 10,text.size = 0,pt.size = 2)
#subset non-B10 breg cells
Breg.nob10=subset(Breg.cut,subset=b10=="Non-B10")
##Figure 5F-feature plot of selected genes in non-B10 Breg cells tsneplot----
myfeatureplot(Breg.nob10,genename = c("Il10","Tgfb1","Ebi3","Il12a"),filename = "non-b10",pt.size = 2,min.cutoff = "q0")

###GSEA plot part####
library(clusterProfiler)
library(org.Mm.eg.db)
library(stringr)
library(ggplot2)

Idents(Breg.cut)="b10"
df.diff.up=FindMarkers(Breg.cut,logfc.threshold = 0,min.pct = 0,ident.1 = "Non-B10")
if(species=="mouse"){
  library(biomaRt)
  human=readRDS("/Users/jalon/Desktop/newest scripts/database/human_gene_ensembl.rds")
  mouse=readRDS("/Users/jalon/Desktop/newest scripts/database/mouse_gene_ensembl.rds")
  total_matrix=df.diff.up
  gene_mouse=row.names(total_matrix)
  print("Converting mouse gene symbols to human...")
  gene_human=getLDS(attributes = c("mgi_symbol"),
                    filters = "mgi_symbol",
                    values=gene_mouse,mart = mouse,
                    attributesL = c("hgnc_symbol"),
                    martL = human,uniqueRows = T)
  gene_human=gene_human[!duplicated(gene_human$MGI.symbol),]
  gene_human=gene_human[!duplicated(gene_human$HGNC.symbol),]
  total_matrix=total_matrix[gene_human$MGI.symbol,]
  row.names(total_matrix)=gene_human$HGNC.symbol
}
gene=rownames(total_matrix)
df=total_matrix
genelist=df$avg_log2FC
names(genelist)=rownames(df)
genelist=sort(genelist,decreasing = T)
#import selected gene sets
gmt=read.gmt("/Users/jalon/Desktop/help_analysis/YSY_B/20210326/B10 VS non-B10.gmt")
gsea=GSEA(geneList = genelist,TERM2GENE = gmt,minGSSize = 6,eps = 0,pvalueCutoff = 1)

##Figure 5G-activated and suppressed GSEA pathways in non-B10 Breg cells----
pdf("GSEA_breg_vs_b10.pdf",width = 1*16,height = 6)
dotplot(gsea,color="pvalue",showCategory = 10,split=".sign")+facet_grid(~.sign)
dev.off()

#construct present Breg markers gene set
Breg_present_marker=c("Cr2","Fcer2a","Cd1d1","Cd5","Il10","Ighd",
                      "Ighm","Sdc1","Cd44","Cd24a","Cxcr4",
                      "Cd22","Ighv7-3","Tnfrsf13b","Havcr1","Il2ra","Cd81","Cd69","Lap3",
                      "Lag3","Igha","Cd274","Faslg","Cd38")
Breg_present_marker=intersect(Breg_present_marker,rownames(Breg.cut))
##Figure S2-feature plot of present Breg cell markers in all B cells tsneplot----
myfeatureplot(all.singlet,genename = Breg_present_marker,filename = "breg_present_marker",pt.size = 0.5,
              min.cutoff = "q0")



