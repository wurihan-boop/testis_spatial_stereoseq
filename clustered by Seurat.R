
#加载空转数据
mamba activate st

rm(list=ls())

library(SeuratDisk)
library(Seurat)
library(ggplot2)
library(viridis)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk",force = TRUE)
getwd()
##1. Data preprocessing
st<- readRDS("/home/wurihan/wrh/data/st_data_rds/DT_A03700A6_bin50_seurat.rds")

plot1 <- VlnPlot(st, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- FeaturePlot(st, features = "nCount_Spatial") + theme(legend.position = "right")

wrap_plots(plot1, plot2)

##normalization
library(future)  
options(future.globals.maxSize = Inf) 
##查看数据结构
> st@
  st@assays        st@meta.data     st@active.assay  st@active.ident  st@graphs        st@neighbors     st@reductions    st@images        st@project.name  st@misc          st@version       st@commands      st@tools
> st@assays$Spatial

st <- SCTransform(st, assay = "Spatial", verbose = FALSE)

p1<-FeaturePlot(st, features = c("DDX4"),reduction="spatial",pt.size=0.5,label=FALSE)
p1
p2<-FeaturePlot(st, features = c("VIM"),reduction="spatial",pt.size=0.5,label=FALSE)
p2

##
st <- RunPCA(st, assay = "SCT", verbose = FALSE)
ElbowPlot(st)

st <- FindNeighbors(st, reduction = "pca", dims = 1:20)
st <- FindClusters(st, resolution = 0.3,verbose = FALSE)
st <- RunUMAP(st, reduction = "pca", dims = 1:20)

table(st@meta.data$seurat_clusters)

p1 <- DimPlot(st, reduction = "umap", label = TRUE)
p2 <- DimPlot(st, label = TRUE,reduction="spatial",label.size = 3)
p1 + p2

save(st,file="/home/wurihan/wrh/st_analysis/DT/DT_A03700A6_bin50_dim20r0.3_clustered.rda")

