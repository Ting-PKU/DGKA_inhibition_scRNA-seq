setwd('~/projects/tumor/')
tmp <- .libPaths()
.libPaths <- c('/home/pengt/projects/tools/anaconda3/envs/R4_2/lib',
               tmp)
library(ggplot2)
library(Seurat)
get_count <- function(countfile) {
  pbmc <- Read10X(data.dir = paste0("data/",countfile))
  pbmc <- CreateSeuratObject(counts = pbmc,  min.cells = 3, min.features = 200)
  pbmc$sampleID <- countfile
  print(countfile)
  return(pbmc)
}
count_ls <- lapply(c('C1','C2','C3','R1','R2'), get_count)
y <- c(count_ls[[2]])
for(i in 3:length(count_ls)){
  y <- c(y, count_ls[[i]])
}
seu <- merge(x = count_ls[[1]],y = y,
             add.cell.ids = c('C1','C2','C3','R1','R2'))
seu <- JoinLayers(seu)
meta <- seu@meta.data
meta$Group <- 'none'
meta[meta$sampleID%in%paste0('C',1:3),]$Group <- 'Control'
meta[meta$sampleID%in%paste0('R',1:2),]$Group <- 'Treatment'
seu$Group <- meta$Group
#### screen cells
VlnPlot(t, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"))
t <- subset(seu, subset = nFeature_RNA > 400
            & nFeature_RNA < 7500 & percent.mt < 30 &nCount_RNA < 50000)
# umap
seu <- NormalizeData(seu)
seu <- FindVariableFeatures(seu, selection.method = "vst",nfeatures = 3000)
seu <- ScaleData(seu, features = rownames(seu))
seu <- RunPCA(seu, features = VariableFeatures(object = seu), npcs = 30)
ElbowPlot(seu)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE)
seu <- RunHarmony(seu, "sampleID", sigma = 0.1, dims.use = 1:20,
                  max.iter.harmony = 15,project.dim = F)
seu <- RunUMAP(seu, dims = 1:20, verbose = FALSE, reduction = "harmony")
DimPlot(seu, group.by = "sampleID", label = F)+NoLegend()
seu <- FindNeighbors(seu, dims = 1:20, verbose = FALSE, reduction = "harmony")
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
DimPlot(seu, group.by = "seurat_clusters", label = TRUE,
        order = T)

markers <- c('GRCm39-Pecam1','GRCm39-Rgs5','GRCm39-Acta2',
             'GRCm39-Epcam','GRCm39-Dcn','GRCm39-Mmp2',
             'GRCm39-Ptprc','GRCh38-EPCAM','GRCh38-PTPRC')
FeaturePlot(seu,'GRCm39-Pecam1')
FeaturePlot(seu,'GRCm39-Epcam')
FeaturePlot(seu,'GRCm39-Dcn')
FeaturePlot(seu,'GRCm39-Mmp2')
FeaturePlot(seu,'GRCm39-Rgs5')
FeaturePlot(seu,'GRCm39-Ptprc')
FeaturePlot(seu,'GRCh38-PTPRC')
FeaturePlot(seu,'GRCh38-PECAM1')
FeaturePlot(seu,'GRCh38-EPCAM')
FeaturePlot(seu,'GRCh38-DCN')
FeaturePlot(seu,'GRCh38-CD14')

seu <- RenameIdents(object = seu, "13" = "Mouse_Endothelial",
                    "14" = "Mouse_Pericyte",
                    "15" = "Mouse_Epithelial", "3" = "Mouse_Fibroblast",
                    "7" = 'Mouse_Fibroblast',"10" = 'Mouse_Fibroblast',
                    "5" = "Mouse_Immune","6" = "Mouse_Immune",
                    "11" = "Mouse_Immune","12" = "Mouse_Immune",
                    "0" = "Human_Epithelial", "1" = "Human_Epithelial",
                    "2" = "Human_Epithelial", "4" = "Human_Epithelial",
                    "9" = "Human_Epithelial",
                    "8" = "Human_Immune","16" = "Human_Immune")
seu$cluster <- Idents(seu)
cols <- c("Human_Epithelial"="#f38989","Mouse_Immune"="#67a4cc","Mouse_Pericyte"="#af93c4",
          "Mouse_Epithelial"="#f9b769","Mouse_Fibroblast"="#B4B549","Human_Immune"="#D38576",
          "Mouse_Endothelial"="#4E856F")
pdf('figure/mix_umap.pdf',width = 3.8,height = 3.8)
DimPlot(seu, order = F,label = T, cols = cols,
        raster = T,pt.size = 1.6)&theme(axis.line = element_blank(),
                                         axis.ticks = element_blank(),
                                         axis.text = element_blank(),
                                         panel.border = element_rect(fill = NA,color = "black",
                                                                     size=1,linetype = "solid"))&
  xlab('UMAP1')&ylab('UMAP2')&NoLegend()
dev.off()
library(wesanderson)
pal <- wes_palette("Zissou1", 50, type = "continuous")
pdf('figure/mix_marker.pdf',width = 6.5,height = 4.1)
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

tmp <- data.frame(seu@reductions$umap@cell.embeddings)
seu$umap1 <- tmp$umap_1
seu$umap2 <- tmp$umap_2
all_seu <- CreateSeuratObject(counts = seu@assays$RNA@layers$counts,  meta.data = seu@meta.data,min.cells = 3, min.features = 200)
saveRDS(all_seu,'source_data/all_seu.rds')
saveRDS(seu@meta.data,'source_data/all_meta.rds')
###### human cells
seu <- subset(seu,subset = cluster%in%c("Human_Immune","Human_Epithelial"))
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
seu <- FindClusters(seu, resolution = 0.3, verbose = FALSE)
DimPlot(seu, group.by = "seurat_clusters", label = TRUE,
        order = T)


t1 <- FindMarkers(seu, ident.1 =c(8),assay = "RNA", min.pct = 0.25, 
                  logfc.threshold = 0.25, test.use = "wilcox",
                  max.cells.per.ident = 2000)
FeaturePlot(seu,'GRCh38-CD3D')
FeaturePlot(seu,'GRCh38-CD14')
seu <- subset(seu,subset = seurat_clusters!=12)
library(RColorBrewer)
pdf('figure/CD3D_featureplot.pdf',width = 4.5,height = 4.2)
FeaturePlot(seu, 'GRCh38-CD3D', raster = T,pt.size = 1.8)&
  scale_colour_gradientn(colours = (brewer.pal(n = 9, name = "Reds")))&
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_rect(fill = NA,color = "black",
                                    size=1.5,linetype = "solid"))&
  xlab('UMAP1')&ylab('UMAP2')
dev.off()
seu <- RenameIdents(object = seu, "8" = "T cells",
                    "9" = "Myeloids",
                    "0" = "Epithelial","1" = "Epithelial",
                    "2" = "Epithelial","3" = "Epithelial",
                    "4" = "Epithelial","5" = "Epithelial",
                    "6" = "Epithelial","7" = "Epithelial",
                    "10" = "Epithelial")
seu$celltype <- Idents(seu)
table(seu$celltype)

saveRDS(seu,'source_data/human_cell_seu.rds')
seu <- readRDS('source_data/human_cell_seu.rds')
seu <- subset(seu,subset=celltype=='T cells')
colors <- c('#D79889',
            '#4E977D','#D5A757','#955E4E')
df <- data.frame(num = as.numeric(table(seu$celltype)),celltype = names(table(seu$celltype)))
pdf('figure/cell_num.pdf',width = 3.6,height = 3.2)
ggplot(df,aes(x = celltype,y = num, fill = celltype,color = celltype))+
geom_bar(stat="identity",width = 0.6)+ 
  xlab('')+ylab('Number of cells')+NoLegend()+
  scale_fill_manual(values = colors)+
  scale_color_manual(values = colors)+
  geom_text(aes(y = num+1000, label = num,), color = "black", size=3.5)+
  mytheme
dev.off()

mytheme <- theme_classic()+theme(legend.position="none",
                                 legend.background = element_rect(fill = NA),
                                 plot.title = element_text(hjust = 0.5),
                                 axis.text.x = element_text(color = 'black'),
                                 axis.text.y = element_text(color = 'black'))
seu$compare <- paste0(seu$Group,'_',seu$celltype)
Idents(seu) <- 'compare'
table(seu$compare)
t1 <- FindMarkers(seu, ident.1 =c('Treatment_T cells'),ident.2 =c('Control_T cells'),assay = "RNA", min.pct = 0.25, 
                  logfc.threshold = 0.25, test.use = "wilcox",
                  max.cells.per.ident = 2000)

EnhancedVolcanoDESeq2 <- function(toptable)
{
  toptable$Significance <- "NS"
  toptable[toptable$p_val<0.01&abs(toptable$avg_log2FC)>0.4,]$Significance <- "sig"
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
    ggtitle('T cells Treatment vs Control')+
    xlab(bquote(~log[2]~"FoldChange")) +
    ylab(bquote(~-log[10]~"p-value")) +
    geom_text_repel(data=subset(toptable, id==T ),
                    aes(label=subset(toptable, id==T )$geneID),
                    size=2,vjust = 0,max.overlaps = 50)+
    geom_vline(xintercept=0.4, 
               linetype="longdash", colour="black", linewidth=0.4)+
    geom_vline(xintercept=-0.4, 
               linetype="longdash", colour="black", linewidth=0.4)+
    geom_hline(yintercept=2, 
               linetype="longdash", colour="black", linewidth=0.4)
  return(plot)
}

library(ggrepel)
t1$gene <- row.names(t1)
t1 <- t1[-grep('GRCm39',t1$gene),]
t1$gene <- gsub('GRCh38-','',t1$gene)
row.names(t1) <- t1$gene
t1$id <- F
t1[c('ISG15','IFIT3','IFI44L','ISG20','IFI44','IFIT1','MX2'),]$id <- T
t1$geneID <- row.names(t1)
pdf('figure/T_vocano.pdf',width = 3.2,height=3.66)
EnhancedVolcanoDESeq2(toptable=t1)
dev.off()

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
t <- go_Enrich(row.names(t1[t1$avg_log2FC<0,]), 'down')
t$LogP <- -log10(t$p.adjust)
t$Description <- factor(t$Description, levels = rev(t$Description))
row.names(t) <-seq(1,dim(t)[1])
df <- t[c(1,4,5,7,9,10,14,15,26),]
p <- ggplot(df,aes(x =Count, y = Description, fill = LogP)) + 
  geom_col(width = 0.85) +theme_bw()
p1 <- p + mytheme+scale_x_continuous(expand = c(0,0))+
  geom_text(data = df,
            aes(x = 0.15, y = Description, label = Description),
            size = 3.8,
            hjust = 0)+
  labs(x = 'Number of genes', y = ' ') + 
  scale_fill_gradientn(colors = (brewer.pal(9,'Spectral')))+
  ggtitle('T cells Treatmnet vs Control up-DEGs')+CenterTitle()+
  theme(legend.position = 'right',axis.text.y = element_blank())
pdf('figure/T_enrich_pathway_up.pdf',width = 5,height=3.66)
p1
dev.off()
