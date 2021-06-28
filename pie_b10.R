pie_b10=function(seurat.object,filename="",organ,label=T){
  Breg_pc=subset(seurat.object,subset=samples==organ)
  Breg_pc$Breg=ifelse(colnames(Breg_pc)%in%b10.cellname,"B10",
                      ifelse(colnames(Breg_pc)%in%breg.cellname,"Breg_nonb10","non-Breg"))
  fm_pc=as.data.frame(table(Breg_pc$Breg))
  colnames(fm_pc)=c("group","count")
  fm_pc$per=round(fm_pc$count/sum(fm_pc$count),4)*100
  pie=ggplot(fm_pc,aes(x="group",y=count,fill=group))+
    geom_bar(color="black",stat = "identity",position = position_stack())+
    coord_polar(theta = "y",start=0)+
    scale_fill_manual(values=c("lightskyblue","hotpink2","mistyrose"),breaks = levels(fm_pc$group))+
    theme_void()
  
  if(label){
    pie=pie+geom_text(aes(y=cumsum(rev(count))-0.5*rev(count),label=paste(rev(per),"%")))
  }
  pdf(paste0("Breg_",filename,"_pie.pdf"), width = 6.5,height = 5.3)
  print(pie)
  dev.off()
}