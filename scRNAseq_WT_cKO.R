library(Seurat)
library(tidyverse)
library(patchwork)
library(ComplexHeatmap)
library(RColorBrewer)
library(Cairo)
install.packages("scCustomize")
library(scCustomize)
library(viridis)

setwd("~/Dropbox (Dropbox @RU)/20230918_JN_singleCell/keepPhgdhInCKO/getClusters_fromMat_keepPhgdhInCKO") #keep
setwd("~/Dropbox (Dropbox @RU)/20230918_JN_singleCell/rmPhgdhInCKO/getClusters_fromMat_rmPhgdhInCKO") #rmv

setwd("~/Dropbox (Dropbox @RU)/20230918_JN_singleCell/keepPhgdhInCKO/getClusters_fromMat_keepPhgdhInCKO")
setwd("~/Dropbox (Dropbox @RU)/20230918_JN_singleCell/rmPhgdhInCKO/getClusters_fromMat_rmPhgdhInCKO")

hfsc.cko.out <- readRDS('seu_SCT_integrate_post_PhgdhFilter_rmMitoClust_mitoRegres_splitInt_20231115.rds')

hfsc.cko.out <- readRDS('seu_SCT_integrate_post_PhgdhFilter_rmMitoClust_mitoRegres_splitInt_wMitoRp_20240514.rds')

DefaultAssay(object = hfsc.cko.out) <- 'SCT'

DefaultAssay(object = hfsc.cko.out) <- 'RNA'

hfsc.cko.out$nCount_RNA

Idents(hfsc.cko.out) <- "orig.ident"
Idents(hfsc.cko.out) <- "seurat_clusters"
Idents(hfsc.cko.out) <- "genotype"
Idents(hfsc.cko.out) <- "time"
Idents(hfsc.cko.out) <- "Phase"
Idents(hfsc.cko.out) <- "geno_ord"
Idents(hfsc.cko.out) <- "time_ord"
Idents(hfsc.cko.out) <- "d5"
Idents(hfsc.cko.out) <- "d3"
Idents(hfsc.cko.out) <- "uw"
hfsc.cko.out$geno_ord <- factor(hfsc.cko.out$genotype, levels = c("WT", "cKO"))
hfsc.cko.out$orig.ident <- factor(hfsc.cko.out$orig.ident, levels = c("UW_WT_1",'UW_cKO_2', 'D3_WT_7', 'D3_cKO_4', 'D5_WT_8', 'D5_cKO_6'))
hfsc.cko.out$d5 <- factor(hfsc.cko.out$orig.ident, levels = c('D5_WT_8', 'D5_cKO_6'))
hfsc.cko.out$d3 <- factor(hfsc.cko.out$orig.ident, levels = c('D3_WT_7', 'D3_cKO_4'))
hfsc.cko.out$uw <- factor(hfsc.cko.out$orig.ident, levels = c('UW_WT_1', 'UW_cKO_2'))

hfsc.cko.out$time.ord <- factor(hfsc.cko.out$time, levels = c('UW', 'D3', 'D5'))
hfsc.cko.out$time.rev <- factor(hfsc.cko.out$time, levels = c('D5', 'D3', 'UW'))

DimPlot(hfsc.cko.out, reduction = "umap", label = T, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white") + NoLegend()

DimPlot(hfsc.cko.out, reduction = "umap", label = T, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white",
        cols = c("#4682B4", "salmon")) + NoLegend()

VlnPlot(hfsc.cko.out, features = c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo'), pt.size = 0, split.by = 'genotype') 

VlnPlot(hfsc.cko.out, assay = 'RNA', log = F, features = c("percent.ribo"), pt.size = 0, split.by = 'geno_ord', cols = c("#4682B4", "salmon")) #, idents = c('7','5'))

VlnPlot(hfsc.cko.out, assay = 'RNA', log = T, features = c("Top2a", 'Krt15', 'Ly6a'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon") , idents = c('1','5','7'))

VlnPlot_scCustom(hfsc.cko.out, features = "percent.ribo", plot_boxplot = T,
                 colors_use = c("#4682B4", "salmon","#4682B4", "salmon","#4682B4", "salmon"),
                 pt.size = 0, log = F) & NoLegend()

VlnPlot_scCustom(hfsc.cko.out, features = "Klf5", plot_boxplot = T,
                 colors_use = c("#4682B4", "salmon","#4682B4", "salmon","#4682B4", "salmon"),
                 pt.size = 0, log = F) & NoLegend()+ theme(axis.title.y=element_blank(),
                                                           axis.title.x=element_blank(),
                                                           axis.text.x=element_blank(),
                                                           axis.ticks.x=element_blank())


VlnPlot(hfsc.cko.out, assay = 'SCT', log = F, features = c('Klf5'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon")) #, idents = c('4','5','7'))

VlnPlot(hfsc.cko.out, assay = 'RNA', log = T, features = c('percent.ribo'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon") , idents = c('4','5','7'))

VlnPlot(hfsc.cko.out, assay = 'SCT', log = F, features = c('Sox9'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon")) + NoLegend() + theme(axis.title.y=element_blank(),
                                                            axis.title.x=element_blank(),
                                                            axis.text.x=element_blank(),
                                                            axis.ticks.x=element_blank())

FeaturePlot_scCustom(hfsc.cko.out, features = c('Sox9'), pt.size = 0.75, 
                     colors_use = viridis_plasma_light_high, 
                     na_color = "lightgray") #+ NoAxes() + NoLegend()

FeaturePlot(hfsc.cko.out, features = c('Klf5', 'Sox9'), pt.size = 0.75, blend = T,
            blend.threshold = .9)

Cluster_Highlight_Plot(seurat_object = hfsc.cko.out, cluster_name = c("0", "1", "5", "9"), highlight_color = c("navy", "forestgreen", "darkred", "darkviolet")) + NoAxes() + NoLegend() #clusters cKO
Cluster_Highlight_Plot(seurat_object = hfsc.cko.out, cluster_name = c("2", "3", "6", "7","11"), highlight_color = c("lightblue", "lightgreen", "lightcoral", "pink", "gold")) + NoAxes() + NoLegend()#clusters cKO

Cluster_Highlight_Plot(seurat_object = hfsc.cko.out, cluster_name = c("1", "7", "5"), highlight_color = c("forestgreen", "pink", "darkred")) #clusters cell cycle
Cluster_Highlight_Plot(seurat_object = hfsc.cko.out, cluster_name = c("1", "7", "5"), highlight_color = c("salmon", "#4682B4", "salmon"))+ NoAxes() + NoLegend() #clusters cell cycle BY GENO

table(hfsc.cko.out$genotype)

isr_genes <- read.csv('isr_genes.csv', header = F)
upr_genes <- read.csv('upr_genes.csv', header = T)
gcn2_genes <- read.csv('response_gcn2_genes.csv', header = T)


#upr_genes
hfsc.cko.out <- AddModuleScore(hfsc.cko.out,  isr_genes, name = 'isr_score')

VlnPlot(hfsc.cko.out, features = 'isr_score1', pt.size = 0, split.by = 'geno_ord') 

VlnPlot(hfsc.cko.out, assay = 'RNA', log = F, features = c('isr_score1'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon")) #, idents = c('4','5','7'))

hfsc.cko.out <- AddModuleScore(hfsc.cko.out,  upr_genes, name = 'upr_score')

VlnPlot(hfsc.cko.out, assay = 'RNA', log = F, features = c('upr_score1'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon")) #, idents = c('4','5','7'))

hfsc.cko.out <- AddModuleScore(hfsc.cko.out,  gcn2_genes, name = 'gcn2_genes')

VlnPlot(hfsc.cko.out, assay = 'RNA', log = F, features = c('gcn2_genes1'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon")) #, idents = c('4','5','7'))

VlnPlot(hfsc.cko.out, features = 'gcn2_genes1', pt.size = 0, split.by = 'genotype') 

VlnPlot(hfsc.cko.out, assay = 'RNA', log = T, features = c('Il24'), pt.size = 0, split.by = 'geno_ord',
        cols = c("#4682B4", "salmon") , idents = c('7','5'))

################# elaine questions 
DimPlot(hfsc.cko.out, reduction = "umap", label = T, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white",
        cols = c("#4682B4", "salmon"), split.by = 'geno_ord') + NoLegend()

DimPlot(hfsc.cko.out, reduction = "umap", label = T, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white",
        cols = c("#4682B4", "salmon"), split.by = 'd5') + NoLegend()

DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white", 
        split.by = 'd5') + NoLegend()

DimPlot(hfsc.cko.out, reduction = "umap", label = T, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "black") + NoLegend() + NoAxes()

DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
        pt.size = 1.25, label.size = 7, label.box = T, label.color = "white",
        split.by = 'geno_ord') + NoLegend() 

DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
        pt.size = 1.25, label.size = 12, label.box = T, label.color = "white",
        cols = c("#4682B4", "salmon")) + NoLegend() + NoAxes()

DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
        pt.size = 1.25, label.size = 12, label.box = T, label.color = "white",
        cols = c("#B79F00", "magenta", '#00BA38')) #+ NoLegend() + NoAxes()

DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
        pt.size = 2, label.size = 12, label.box = T, label.color = "white",
        cols = c("#4682B4", "salmon"), split.by = 'time.rev') + NoLegend() 
table(Idents(hfsc.cko.out))


leiden_plot <- DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
                       pt.size = 1.25, label.size = 7, label.box = T, label.color = "black") + NoLegend() + NoAxes()

leiden_plot
DimPlot(hfsc.cko.out, reduction = "umap", label = F, repel = T, 
                       pt.size = 1.25, label.size = 7, label.box = T, label.color = "black") + NoAxes()

LabelClusters(leiden_plot, id = "ident",  size = 8, box = T, fontface = "bold", color = "white")


