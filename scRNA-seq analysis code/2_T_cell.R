setwd('~/projects/tumor/')
library(ggplot2)
library(Seurat)
library(harmony)
DimPlot(seu)
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures = 3000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu), npcs = 30)
ElbowPlot(seu)
seu <- RunHarmony(seu, "sampleID", sigma = 0.1, dims.use = 1:15,
                  max.iter.harmony = 15,project.dim = F)
seu <- RunUMAP(seu, dims = 1:15, verbose = FALSE, reduction = "harmony")
DimPlot(seu, group.by = "sampleID", label = F)+NoLegend()
seu <- FindNeighbors(seu, dims = 1:15, verbose = FALSE, reduction = "harmony")
seu <- FindClusters(seu, resolution = 1, verbose = FALSE)
DimPlot(seu, group.by = "seurat_clusters", label = TRUE,
        order = T)

t1 <- FindMarkers(seu, ident.1 =c(3),assay = "RNA", min.pct = 0.25, 
                  logfc.threshold = 0.25, test.use = "wilcox",
                  max.cells.per.ident = 2000)
FeaturePlot(seu,'GRCh38-CD3D')
FeaturePlot(seu,'GRCh38-CD8A')
FeaturePlot(seu,'GRCh38-CD4')
FeaturePlot(seu,'GRCh38-FOXP3')
FeaturePlot(seu,'GRCh38-MKI67')
FeaturePlot(seu,'GRCh38-TIGIT')
FeaturePlot(seu,'GRCh38-CD4')
FeaturePlot(seu,'GRCh38-NCAM1')

FeaturePlot(seu,'GRCh38-TCF7')
FeaturePlot(seu,'GRCh38-SELL')
FeaturePlot(seu,'GRCh38-IL7R')
FeaturePlot(seu,'GRCh38-GZMH')
seu <- subset(seu,subset = seurat_clusters!=12)
library(RColorBrewer)

seu <- RenameIdents(object = seu, "5" = "MKI67+ T cells",
                    "7" = "Tnaive",
                    "1" = "Tem","4" = "Tem",
                    "0" = "Tem","2" = "Teff",
                    "6" = "Tex","3" = "Treg","8" = "NK")
seu$celltype <- Idents(seu)
DimPlot(seu,group.by = 'celltype',label = T)
table(seu$celltype)

saveRDS(seu,'source_data/human_T_cell_seu.rds')

seu <- readRDS('source_data/human_T_cell_seu.rds')

cols <- c("MKI67+ T cells"="#67a4cc","Tnaive"="#B4B549","Tem"="#af93c4",
          "Teff"="#f38989","Tex"="#f9b769","Treg"="#D38576",
          "NK"="#4E856F")
pdf('figure/T_umap.pdf',width = 3.8,height = 3.8)
DimPlot(seu, order = F,label = T, cols = cols,
        raster = T,pt.size = 5)&theme(axis.line = element_blank(),
                                        axis.ticks = element_blank(),
                                        axis.text = element_blank(),
                                        panel.border = element_rect(fill = NA,color = "black",
                                                                    size=1,linetype = "solid"))&
  xlab('UMAP1')&ylab('UMAP2')&NoLegend()
dev.off()


library(wesanderson)
pal <- wes_palette("Zissou1", 50, type = "continuous")
pdf('figure/T_cell_marker.pdf',width = 6,height = 3.8)
markers <- c('MKI67','TCF7','IL7R','GZMH','GZMK',
             'TIGIT','FOXP3','NCAM1')
markers <- paste0('GRCh38-',markers)
DotPlot(seu, features = markers, dot.scale = 10) + RotatedAxis() +
  theme_classic()+theme(legend.position="right",
                        legend.background = element_rect(fill = NA),
                        axis.line = element_line(size = 0),
                        panel.border = element_rect(size = 1,fill=NA,color = 'black'),
                        plot.title = element_text(hjust = 0.5),
                        axis.text.x = element_text(angle = 30, hjust=1,color = 'black'),
                        axis.text.y = element_text(color = 'black',face = 'bold'))+
  scale_colour_gradientn(colours = pal)+ 
  labs(y = "", x="")
dev.off()


pdf('figure/T_umap_split.pdf',width = 6.5,height = 3.8)
DimPlot(seu, order = F,label = T,split.by = 'Group', cols = cols,
        raster = T,pt.size = 5)&theme(axis.line = element_blank(),
                                      axis.ticks = element_blank(),
                                      axis.text = element_blank(),
                                      panel.border = element_rect(fill = NA,color = "black",
                                                                  size=1,linetype = "solid"))&
  xlab('UMAP1')&ylab('UMAP2')&NoLegend()
dev.off()

meta <- seu@meta.data
meta$cellid <- row.names(meta)
t <- aggregate(meta$cellid,by=list(sample = meta$sampleID),length)
tt <- aggregate(meta$cellid,by=list(sample = meta$sampleID,
                                    celltype = meta$celltype),length)

tmp <- merge(t,tt,by='sample')
tmp$ratio <- tmp$x.y/tmp$x.x
tmp$Group <- 'Control'
tmp[tmp$sample%in%c('R1','R2'),]$Group <- 'Treatment'

tt <- aggregate(tmp$ratio,by=list(Group = tmp$Group,
                                celltype = tmp$celltype),mean)
t <- aggregate(tt$x,by=list(Group = tt$Group),sum)
tmp <- merge(t,tt,by='Group')
tmp$ratio <- tmp$x.y/tmp$x.x
pdf('figure/T_cluster_freq_barplot.pdf',width = 5.4,height = 4.45)
ggplot(tmp, aes(fill=celltype, y=ratio, x=Group)) +
  geom_bar( position="stack", stat="identity",width = 0.6,
            alpha = 1) + 
  scale_fill_manual(values = cols)+
  theme_bw()+theme(axis.text.x = element_text(angle=30, vjust=0.85, hjust = 0.7))+
  scale_y_continuous(expand=c(0,0)) +
  #geom_text(position = position_stack(vjust = 0.5),
  #         aes(label = data), family = 'Arial') +
  ylab('Mean frequency of clusters') + xlab("")+
  guides(fill=guide_legend(title=""))
dev.off()


seu$compare <- paste0(seu$Group,'_',seu$celltype)
Idents(seu) <- 'compare'
table(seu$compare)
t1 <- FindMarkers(seu, ident.1 =c('Treatment_Teff'),ident.2 =c('Control_Teff'),assay = "RNA", min.pct = 0.25, 
                  logfc.threshold = 0.25, test.use = "wilcox",
                  max.cells.per.ident = 2000)
t1 <- t1[t1$p_val<0.05,]
tt <- t1[-grep('GRCm',row.names(t1)),]
row.names(tt) <- gsub('GRCh38-','',row.names(tt))
tt$Group <- 'Treatment_up'
tt[tt$avg_log2FC<0,]$Group <- 'Treatment_down'
write.csv(tt,'source_data/Tex_diff_genes.csv',quote = F)
EnhancedVolcanoDESeq2 <- function(toptable)
{
  toptable$Significance <- "NS"
  toptable[toptable$p_val<0.05&abs(toptable$avg_log2FC)>0.4,]$Significance <- "sig"
  print(table(toptable$Significance))
  toptable$Significance <- factor(toptable$Significance, levels=c("NS", "sig"))
  plot <- ggplot(toptable, aes(x=avg_log2FC, y=-log10(p_val))) +
    geom_point(aes(color=factor(Significance)), alpha=1/2, size=1.8) +
    scale_color_manual(values=c(NS="grey30",  sig="#F93838")) +
    mytheme+theme(legend.position = "none",
                  axis.line = element_line(size = 0),
                  panel.border = element_rect(linewidth = 1,
                                              fill=NA,color = 'Black'),
                  plot.title = element_text(hjust = 0.5))+
    xlab(bquote(~log[2]~"FoldChange")) +
    ylab(bquote(~-log[10]~"p-value")) +
    geom_text_repel(data=subset(toptable, id==T ),
                    aes(label=subset(toptable, id==T )$geneID),
                    size=3,vjust = 0,max.overlaps = 50)+
    geom_vline(xintercept=0.4, 
               linetype="longdash", colour="black", linewidth=0.4)+
    geom_vline(xintercept=-0.4, 
               linetype="longdash", colour="black", linewidth=0.4)+
    geom_hline(yintercept=-log10(0.05), 
               linetype="longdash", colour="black", linewidth=0.4)
  return(plot)
}

library(ggrepel)
t1 <- FindMarkers(seu, ident.1 =c('Treatment_Teff'),ident.2 =c('Control_Teff'),assay = "RNA", min.pct = 0.25, 
                  logfc.threshold = 0.25, test.use = "wilcox",
                  max.cells.per.ident = 2000)
t1$gene <- row.names(t1)
t1 <- t1[-grep('GRCm39',t1$gene),]
t1$gene <- gsub('GRCh38-','',t1$gene)
t1 <- t1[-grep('ENSG',t1$gene),]
row.names(t1) <- t1$gene
t1$id <- F
t1[c('OAS1','OAS3','ISG15','ISG20','ICOS','THEMIS','MX2','IFI35'),]$id <- T
t1$geneID <- row.names(t1)
pdf('figure/Teff_vocano.pdf',width = 3.2,height=3.66)
EnhancedVolcanoDESeq2(toptable=na.omit(t1))+
  ggtitle('Teff cells Treatment vs Control')
dev.off()
meta <- seu@meta.data
meta$gene <- seu@assays$RNA$data['GRCh38-TIGIT',]
meta <- meta[meta$celltype=='Tex',]
library(rstatix,lib.loc = '~/projects/tools/anaconda3/envs/R4_2/lib/R/library/')
library(ggpubr,lib.loc = '~/projects/tools/anaconda3/envs/R4_2/lib/R/library/')
stat.test <- meta %>%
  wilcox_test(gene ~ Group) %>%
  add_significance()%>% add_xy_position(x = c("Group"))
stat.test
data <- meta
p <- ggplot(data,aes(y=gene,x=Group))+
  geom_violin( aes(fill=Group),color = 'white', width=0.8, size=0.5, scale="width")+
  geom_boxplot(size = 0.05,width = 0.05,aes(colour = Group),outlier.shape = NA) +
  theme_pubr() +
  theme(axis.title.y=element_text(size=12),
        axis.text.y =element_text(size = 10)) + 
  scale_fill_manual(values = c("#8cd17d","#f1cd63")) +
  scale_colour_manual(values = c("#8cd17d","#f1cd63"))+NoLegend()+
  xlab('')+ylab('Expression level of TIGIT')+
  stat_pvalue_manual(stat.test, label = "p.signif",
                     bracket.size = 0.75,size = 8)
pdf('figure/TIGIT_compare_vocano.pdf',width = 3.7,height = 3.9)
p
dev.off()
FeaturePlot(seu, 'GRCh38-ENTPD1',split.by = 'Group')&
  scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))&
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(fill = NA,color = "black",
                                    size=1.5,linetype = "solid"))&
  xlab('UMAP1')&ylab('UMAP2')

library(clusterProfiler)
go_Enrich <- function(Gene_list,Group){
  tmp <- enrichGO(gene = Gene_list,
                  OrgDb = 'org.Hs.eg.db',
                  keyType = 'SYMBOL',
                  ont = 'BP',
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1)
  tmp <- data.frame(tmp)[,c(2,6,8,9)]
  tmp$Group <- Group
  return(tmp)
}
t <- go_Enrich(row.names(t1[t1$avg_log2FC>0,]), 'up')
t$LogP <- -log10(t$p.adjust)
t$Description <- factor(t$Description, levels = rev(t$Description))
row.names(t) <-seq(1,dim(t)[1])
df <- t[c(1,2,3,4,7,10,21,27,28),]
p <- ggplot(df,aes(x =Count, y = Description, fill = LogP)) + 
  geom_col(width = 0.85) +theme_bw()
p1 <- p + mytheme+scale_x_continuous(expand = c(0,0))+
  geom_text(data = df,
            aes(x = 0.15, y = Description, label = Description),
            size = 3.8,
            hjust = 0)+
  labs(x = 'Number of genes', y = ' ') + 
  scale_fill_gradientn(colors = rev(brewer.pal(9,'Spectral')))+
  ggtitle('Teff cells Treatmnet vs Control up-DEGs')+CenterTitle()+
  theme(legend.position = 'right',axis.text.y = element_blank())
pdf('figure/Teff_enrich_pathway_up.pdf',width = 5.1,height=3.66)
p1
dev.off()
FeaturePlot(seu,'GRCh38-TIGIT',split.by = 'Group')
