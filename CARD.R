########## CARD 空间转录组去卷积分析 ##########################
# Author: 
# Date: 
# Description: CARD deconvolution for donkey ST bin80 data (DT), reference: snRNA-seq (snRNA.rda, Seurat object named pbmc)
# Requirement: Seurat, CARD, MuSiC, RColorBrewer, ggplot2, SingleCellExperiment
# Usage: Modify paths in the first section before running
##############################################################

## ========== 用户修改区（只改这里） ==========
set.seed(12345)  # 固定随机种子，保证去卷积可复现
root_dir     <- "/home/wurihan/wrh/st_analysis/DT/DT_bin80_CARD/"
st_rda_path  <- "/home/wurihan/wrh/st_analysis/DT/DT_A03700A6_bin80_dim20r0.5_clustered.rda"
sc_rda_path  <- "snRNA.rda"
ncore_use    <- 10  # CARD_SCMapping并行核数，不要直接40，根据服务器调整
## ===========================================

# 加载包
library(Seurat)
library(CARD)
library(MuSiC)
library(RColorBrewer)
library(ggplot2)
library(SingleCellExperiment)

# 创建输出目录
dir.create(root_dir, recursive = TRUE, showWarnings = FALSE)
setwd(root_dir)
data_dir    <- file.path(root_dir, "DT_ref_donkey_bin80_CARD_data/")
figure_dir  <- file.path(root_dir, "DT_ref_donkey_bin80_CARD_figure/")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
cat("Working directory: ", getwd(), "\n")

##### 1. Create CARD object
## ST data 处理
load(st_rda_path)
spatial_count <- st@assays$SCT@counts
spatial_location <- st@reductions$spatial@cell.embeddings
spatial_location <- as.data.frame(spatial_location)
colnames(spatial_location) <- c("x", "y")

## snRNA reference 单细胞参考数据（对象名 pbmc）
load(sc_rda_path)
sc_count  <- pbmc@assays$RNA@counts
sc_meta   <- pbmc@meta.data

# 全局统一细胞配色
celltype_colors <- c(
  "SPG" = "#7370dd",
  "SPC" = "#ff946b",
  "Spermatid" = "#00b8f6",
  "Sertoli cell" = "#cf909c",
  "Leydig cell" = "#4169E1",
  "Myoid cell" = "#88c992",
  "Endothelial cell" = "#CD5C5C",
  "Immune cell" = "#808080"
)

### 构建CARD对象
CARD_obj <- createCARDObject(
  sc_count = sc_count,
  sc_meta = sc_meta,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  ct.varname = "celltype",
  ct.select = unique(sc_meta$celltype),
  sample.varname = "orig.ident",
  minCountGene = 100,
  minCountSpot = 5
)

########## 2. Deconvolution using CARD
CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
print(head(CARD_obj@Proportion_CARD, 2))
save(CARD_obj, file = file.path(data_dir, "DT_CARD_bin80_deconvolution_result.rda"))

######### 3. 可视化 ########################################
## 3.1 饼图：每个spot细胞类型比例
p1 <- CARD.visualize.pie(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  colors = celltype_colors,
  radius = NULL
)
print(p1)
ggsave(file.path(figure_dir, "DT_CARD_bin80_Proportion_CARD_pie.pdf"),
       plot = p1, width = 20, height = 20, bg = "white")

## 3.2 单类型示例：SPC
ct.visualize <- c("SPC")
p2_spc <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  ct.visualize = ct.visualize,
  colors = c("white", celltype_colors["SPC"]),
  NumCols = 1,
  pointSize = 0.5
)
print(p2_spc)
ggsave(file.path(figure_dir, "DT_CARD_bin80_distribution_of_SPC.pdf"),
       plot = p2_spc, width = 10, height = 10, bg = "white")

## 3.3 单类型示例：SPG
ct.visualize <- c("SPG")
p2_spg <- CARD.visualize.prop(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  ct.visualize = ct.visualize,
  colors = c("white", celltype_colors["SPG"]),
  NumCols = 1,
  pointSize = 0.5
)
print(p2_spg)
ggsave(file.path(figure_dir, "DT_CARD_bin80_distribution_of_SPG.pdf"),
       plot = p2_spg, width = 10, height = 10, bg = "white")

## 3.4 循环批量绘制所有细胞类型（使用自定义配色）
cell_types <- colnames(CARD_obj@Proportion_CARD)
for (ct in cell_types) {
  p_ct <- CARD.visualize.prop(
    proportion = CARD_obj@Proportion_CARD,
    spatial_location = CARD_obj@spatial_location,
    ct.visualize = ct,
    colors = c("white", celltype_colors[ct]),
    NumCols = 1,
    pointSize = 0.5
  )
  print(p_ct)
  fname <- paste0("DT_bin80_clustered_spatial_distribution_of_", gsub(" ", "_", ct), ".pdf")
  ggsave(file.path(figure_dir, fname), plot = p_ct, width = 10, height = 10, bg = "white")
}

## 3.5 两种细胞共可视化 SPC & SPG
p3 <- CARD.visualize.prop.2CT(
  proportion = CARD_obj@Proportion_CARD,
  spatial_location = CARD_obj@spatial_location,
  ct2.visualize = c("SPC","SPG"),
  colors = list(c("white", celltype_colors["SPC"]), c("white", celltype_colors["SPG"]))
)
print(p3)
ggsave(file.path(figure_dir, "DT_bin80_2C_SPC_SPG.pdf"), plot = p3, width = 20, height = 20, bg = "white")

## 3.6 细胞类型比例相关性热图
p4 <- CARD.visualize.Cor(CARD_obj@Proportion_CARD, colors = NULL)
print(p4)
ggsave(file.path(figure_dir, "DT_bin80_Cor.pdf"), plot = p4, width = 10, height = 10, bg = "white")

############### Refined spatial map ###############
#### 4.1 Imputation on the newly grided spatial locations
CARD_obj <- CARD.imputation(CARD_obj, NumGrids = 2000, ineibor = 10, exclude = NULL)
location_imputation <- cbind.data.frame(
  x = as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop), split="x"),"[",1)),
  y = as.numeric(sapply(strsplit(rownames(CARD_obj@refined_prop), split="x"),"[",2))
)
rownames(location_imputation) <- rownames(CARD_obj@refined_prop)

p5 <- ggplot(location_imputation, aes(x = x, y = y)) +
  geom_point(shape=22, color = "#7dc7f5") +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position="bottom",
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5))
print(p5)
ggsave(file.path(figure_dir, "DT_bin80_location_imputation.pdf"), plot = p5, width = 20, height = 20, bg = "white")

## 4.2 高分辨率插值后，一次性绘制全部细胞类型（修复ct.visualize悬空bug）
color_palette <- c("#00008B", "#1E90FF", "#00FFFF", "#7FFFAA", "#FFFF00", "#FF4500", "#DC143C")
p6 <- CARD.visualize.prop(
  proportion = CARD_obj@refined_prop,
  spatial_location = location_imputation,
  ct.visualize = cell_types,
  pointSize = 0.5,
  colors = color_palette,
  NumCols = 8
) +
  theme_bw() +
  theme(
    panel.background = element_blank(),
    plot.background = element_blank(),
    strip.background = element_blank(),
    panel.border = element_blank(),
    strip.text = element_text(size = 16, face = "bold"),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.box.margin = margin(0, 5, 0, 5)
  )
print(p6)
ggsave(file.path(figure_dir, "DT_bin80_refined_allcelltype.pdf"),
       plot = p6, width=22, height=3, bg="white", device = "pdf", dpi=600)
ggsave(file.path(figure_dir, "DT_bin80_refined_allcelltype.tiff"),
       plot = p6, bg="white", dpi = 600, width=22, height=3, device = "tiff")

##4.3 标记基因高分辨率插值图
# 候选标记基因列表
your_genes <- c(
  "PRM2","HOOK1","TNP1","TPPP2","SPATA3","PRM3",
  "AKAP4","TEX29","TNP2","PRM1","ODF1","ODF2","OAZ3"
)
# 自动过滤refined表达矩阵中不存在的基因
gene_avail <- your_genes[your_genes %in% rownames(CARD_obj@refined_expression)]
cat("Available marker genes for refined plot: ", gene_avail, "\n")

p7 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@refined_expression,
  spatial_location = location_imputation,
  gene.visualize = c("FOXO1","SPDYA","ODF2","SOX9","STAR","DCN","PECAM1","PTPRC"),
  colors = NULL,
  NumCols = 8
)
print(p7)
ggsave(file.path(figure_dir, "DT_bin80_marker_gene_expression_enhanced_resolution.pdf"),
       plot = p7, width=18, bg="white", dpi=300, height=3)
ggsave(file.path(figure_dir, "DT_bin80_marker_gene_expression_enhanced_resolution.tiff"),
       plot = p7, bg="white", dpi=300, width=18, height=3)

##4.4 原始spot水平标记基因图（不插值）
p8_syce3 <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("SYCE3"),
  NumCols = 1
)
print(p8_syce3)
ggsave(file.path(figure_dir, "DT_bin80_SYCE3_spot_level.pdf"),
       plot = p8_syce3, width=10, height=10, bg="white")

p8_8gene <- CARD.visualize.gene(
  spatial_expression = CARD_obj@spatial_countMat,
  spatial_location = CARD_obj@spatial_location,
  gene.visualize = c("FOXO1","SPDYA","ODF2","SOX9","STAR","DCN","PECAM1","PTPRC"),
  colors = NULL,
  NumCols = 8
)
print(p8_8gene)
ggsave(file.path(figure_dir, "DT_bin80_marker_gene_expression_countMat.pdf"),
       plot = p8_8gene, width=18, bg="white", dpi=300, height=3)
ggsave(file.path(figure_dir, "DT_bin80_marker_gene_expression_countMat.tiff"),
       plot = p8_8gene, bg="white", dpi=300, width=18, height=3)

###5. Extension of CARD for single cell resolution mapping
scMapping <- CARD_SCMapping(CARD_obj, shapeSpot="Square", numCell=20, ncore = ncore_use)
print(scMapping)

MapCellCords <- as.data.frame(colData(scMapping))
count_SC <- assays(scMapping)$counts

colors_sc <- c("#8DD3C7","#CFECBB","#F4F4B9","#CFCCCF","#D1A7B9","#E9D3DE","#F4867C","#C0979F",
               "#D5CFD6","#86B1CD","#CEB28B","#EDBC63","#C59CC5","#C09CBF","#C2D567","#C9DAC3","#E1EBA0",
               "#FFED6F","#CDD796","#F8CDDE")
p10 <- ggplot(MapCellCords, aes(x = x, y = y, colour = CT)) +
  geom_point(size = 3.0) +
  scale_colour_manual(values = colors_sc) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_rect(colour = "white", fill="white"),
        plot.background = element_rect(colour = "white", fill="white"),
        legend.position="bottom",
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 13,face="bold"),
        legend.text=element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 15,face="bold"))+
  guides(color=guide_legend(title="Cell Type"))
print(p10)
ggsave(file.path(figure_dir, "DT_bin80_scresolution_allcelltype.pdf"),
       plot = p10, width=20, height=20, bg="white")

