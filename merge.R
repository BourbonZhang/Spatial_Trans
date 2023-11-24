library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)


MS1 <- Seurat::Load10X_Spatial("/media/zwh/DATADRIVE0/work/8001-MS/Sample_8001-MS-1-GEX_CAGGACTCTA-CCATCCTACG")
MS2 <- Seurat::Load10X_Spatial("/media/zwh/DATADRIVE0/work/8001-MS/Sample_8001-MS-2-GEX_CTACTGAATT-AGCGCTAGCG")

Idents(MS1) <- "WT"
Idents(MS2) <- "KO"

# 如果需要，更新细胞名称以确保它们是唯一的
colnames(MS1) <- paste("WT", colnames(MS1), sep = "_")
colnames(MS2) <- paste("KO", colnames(MS2), sep = "_")

# 合并数据
brain.merge <- merge(x = MS1, y = MS2, project = "MS_Comparison")
brain.merge[["orig.ident"]] <- Idents(brain.merge)

# 执行 SCTransform 正则化
brain.merge <- SCTransform(brain.merge, assay = "Spatial", verbose = TRUE)

backup1 <- brain.merge

# 降维、聚类、运行UMAP
brain.merge <- RunPCA(brain.merge, features = VariableFeatures(object = brain.merge))
ElbowPlot(brain.merge)
brain.merge <- FindNeighbors(brain.merge, dims = 1:14)
brain.merge <- FindClusters(brain.merge, resolution = 3.75) # 选择适当的分辨率
brain.merge <- RunUMAP(brain.merge, dims = 1:14)

DimPlot(brain.merge, label = T, reduction = "umap", group.by = c("ident", "orig.ident"), )

SpatialDimPlot(brain.merge, label = F, label.size = 2, pt.size = 1.2, alpha = 1)

#分开的umap
# 分别为WT和KO细胞创建子集Seurat对象
wt_cells <- subset(brain.merge, orig.ident == "WT")
ko_cells <- subset(brain.merge, orig.ident == "KO")

# 为WT细胞绘制UMAP图
DimPlot(wt_cells, reduction = "umap", label = TRUE) + ggtitle("WT Cells")

# 为KO细胞绘制UMAP图
DimPlot(ko_cells, reduction = "umap", label = TRUE) + ggtitle("KO Cells")
#单独的cluster分析
WT_markers <- FindAllMarkers(wt_cells, only.pos = F, min.pct = 0.5, logfc.threshold = 0.25, test.use = "wilcox", pvalue.cutoff = 0.01, core = 4, t = 10)

WT_top10_markers <- as.data.frame(WT_markers %>% group_by(cluster) %>% 
                                top_n(n = 10, wt = avg_log2FC))
WT_top10_markers


KO_markers <- FindAllMarkers(ko_cells, only.pos = F, min.pct = 0.5, logfc.threshold = 0.25, test.use = "wilcox", pvalue.cutoff = 0.01, core = 4)

KO_top10_markers <- as.data.frame(KO_markers %>% group_by(cluster) %>% 
                                    top_n(n = 10, wt = avg_log2FC))
WT_top10_markers

