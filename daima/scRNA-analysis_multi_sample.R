
#### 1. 加载分析使用的工具包 ####
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggsci)
cluster_cols <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                  "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                  "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                  "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                  "#AA8282", "#D4B7B7", "#8600BF", "#BA5CE3", "#808000",
                  "#AEAE5C", "#1E90FF", "#00BFFF", "#56FF0D", "#FFFF00",
                  "#FB9A99", "#E31A1C", "#6A3D9A", "#B15928", "#FDBF6F")
##### 2. 读入原始表达数据 ####
setwd("D:/0502/0526_DKD/GSE209781/GSE209781_RAW")
#10X 数据
DKD01 <- Read10X("./DKD01/")
DKD02 <- Read10X("./DKD02/")
DKD03 <- Read10X("./DKD03/")
NM01 <- Read10X("./NM01/")
NM02 <- Read10X("./NM02/")
NM03 <- Read10X("./NM03/")
# 这里的列名就是barcode
colnames(DKD01) <- paste(colnames(DKD01),"DKD01",sep = "_")
colnames(DKD02) <- paste(colnames(DKD02),"DKD02",sep = "_")
colnames(DKD03) <- paste(colnames(DKD03),"DKD03",sep = "_")
colnames(NM01) <- paste(colnames(NM01),"NM01",sep = "_")
colnames(NM02) <- paste(colnames(NM02),"NM02",sep = "_")
colnames(NM03) <- paste(colnames(NM03),"NM03",sep = "_") 
#将所有读入的数据合并成一个大的矩阵
#合并时需注意行名一致
#既有10X的数据又有表达矩阵的数据，全部转换为表达矩阵再进行合并
#关于矩阵合并请见单独的矩阵合并脚本“merge_matrix.R”
experiment.data <- cbind(DKD01,DKD02,DKD03,NM01,NM02,NM03)
#save(experiment.data,file = 'experiment.data_dgCMatrix.RData')
#创建一个文件夹用于写分析结果
sam.name <- "multi"
if(!dir.exists(sam.name)){
  dir.create(sam.name)
}

#### 3. 创建Seurat分析对象 ####
experiment.aggregate <- CreateSeuratObject(
  experiment.data,
  project = "multi", 
  min.cells = 10,
  min.features = 200,
  names.field = 2,
  names.delim = "-")
save(experiment.aggregate,file = 'experiment.aggregate_before_filter.RData')
#### 4. 数据概览 & QC ####
#查看SeuratObject中的对象
#rm(experiment.aggregate1)
gc()
slotNames(experiment.aggregate)
#assay
experiment.aggregate@assays
#细胞及细胞中基因与RNA数量
dim(experiment.aggregate@meta.data)
View(experiment.aggregate@meta.data)
table(experiment.aggregate$orig.ident)#查看各组细胞数

##QC：统计线粒体基因在每个细胞中的占比
experiment.aggregate[["percent.MT"]] <- PercentageFeatureSet(experiment.aggregate, 
                                                             pattern = "^MT")
pdf(paste0("./",sam.name,"/QC-VlnPlot.pdf"),width = 15,height = 5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), cols = cluster_cols,
        ncol = 3)
dev.off()

##QC：统计基因数，RNA，线粒体基因分布
gene.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nFeature_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
rna.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$nCount_RNA,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
MT.freq <- do.call("cbind", tapply(experiment.aggregate@meta.data$percent.MT,experiment.aggregate@meta.data$orig.ident,quantile,probs=seq(0,1,0.05)))
freq.combine <- as.data.frame(cbind(gene.freq,rna.freq,MT.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq),"Gene",sep = "_"),
                            paste(colnames(rna.freq),"RNA",sep = "_"),
                            paste(colnames(MT.freq),"MT",sep = "_"))
write.table(freq.combine,file = paste0(sam.name,"/QC-gene_frequency.txt"),quote = F,sep = "\t")
rm(gene.freq,rna.freq,MT.freq)
View(freq.combine)

##QC：基因数与线粒体基因以及RNA数量的分布相关性
plot1 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "percent.MT",cols = cluster_cols,)
plot2 <- FeatureScatter(experiment.aggregate, feature1 = "nCount_RNA", feature2 = "nFeature_RNA",cols = cluster_cols,)
pdf(paste0("./",sam.name,"/QC-FeatureScatter.pdf"),width = 8,height = 4.5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()
rm(plot1,plot2)
save(experiment.aggregate,file = 'experiment.aggregate_before_filter1.RData')

#### 5. 筛选细胞 ####
cat("Before filter :",nrow(experiment.aggregate@meta.data),"cells\n")
experiment.aggregate <- subset(experiment.aggregate, 
                               subset = 
                                 nFeature_RNA > 500 & 
                                 nFeature_RNA < 5000 & 
                                 nCount_RNA > 100 & 
                                 nCount_RNA < 30000 &
                                 percent.MT < 30)
cat("After filter :",nrow(experiment.aggregate@meta.data),"cells\n")
table(experiment.aggregate$orig.ident)#查看过滤后各组细胞数
##过滤后的质控图
pdf(paste0("./",sam.name,"/QC-VlnPlot_after.pdf"),width = 15,height =5)
VlnPlot(experiment.aggregate, features = c("nFeature_RNA", "nCount_RNA", "percent.MT"), cols = cluster_cols,
        ncol = 3)
dev.off()

#### 6. 表达量标准化 ####
experiment.aggregate <- NormalizeData(experiment.aggregate, 
                                      normalization.method = "LogNormalize",
                                      scale.factor = 10000)

#计算表达量变化显著的基因FindVariableFeatures
experiment.aggregate <- FindVariableFeatures(experiment.aggregate, 
                                             selection.method = "vst",
                                             nfeatures = 2000)

#展示标准化之后的整体表达水平
top10 <- head(x = VariableFeatures(experiment.aggregate), 10)
plot1 <- VariableFeaturePlot(experiment.aggregate)
plot2 <- LabelPoints(plot = plot1, points = top10)
pdf(file = paste0(sam.name,"/Norm-feature_variable_plot.pdf"),width = 8,height = 5)
CombinePlots(plots = list(plot1, plot2),legend = "none")
dev.off()

#### 7. 均一化与PCA ####
#均一化（需要一点时间）
experiment.aggregate <- ScaleData(
  object = experiment.aggregate,
  do.scale = FALSE,
  do.center = FALSE,
  vars.to.regress = c("orig.ident","percent.MT"))

#PCA计算
experiment.aggregate <- RunPCA(object = experiment.aggregate, 
                               features = VariableFeatures(experiment.aggregate),
                               verbose = F,npcs = 50)

#PCA结果展示-1
pdf(paste0("./",sam.name,"/PCA-VizDimLoadings.pdf"),width = 7,height = 7)
VizDimLoadings(experiment.aggregate, dims = 1:2, reduction = "pca")
dev.off()

#PCA结果展示-2
pdf(paste0("./",sam.name,"/PCA-DimPlot.pdf"),width = 5,height = 4)
DimPlot(experiment.aggregate, reduction = "pca",cols = cluster_cols,)
dev.off()

#PCA结果展示-3
pdf(paste0("./",sam.name,"/PCA-DimHeatmap.pdf"),width = 5,height = 4)
DimHeatmap(experiment.aggregate, dims = 1:6, cells = 500, balanced = TRUE)
dev.off()

#### 8. 确定细胞类群分析PC ####
#碎石图
pdf(paste0("./",sam.name,"/PCA-ElbowPlot.pdf"),width = 6,height = 5)
ElbowPlot(experiment.aggregate,ndims = 40)
dev.off()

#确定用于细胞分群的PC
dim.use <- 1:30

#### 9. 细胞分群TSNE算法 ####
#TSNE算法
experiment.aggregate <- FindNeighbors(experiment.aggregate, dims = dim.use)
experiment.aggregate <- FindClusters(experiment.aggregate, resolution = 0.5)
experiment.aggregate <- RunTSNE(experiment.aggregate, dims = dim.use, 
                                do.fast = TRUE)
pdf(paste0("./",sam.name,"/CellCluster-TSNEPlot_res0.5_",max(dim.use),"PC.pdf"),width = 6,height = 5)
DimPlot(object = experiment.aggregate, pt.size=0.5,label = T,cols = cluster_cols)
dev.off()

table(experiment.aggregate@meta.data$orig.ident)

##### 10. 计算marker基因 #####
#这一步计算的时候可以把min.pct以及logfc.threshold调的比较低，然后再基于结果手动筛选
all.markers <- FindAllMarkers(experiment.aggregate, only.pos = TRUE, 
                              min.pct = 0.3, logfc.threshold = 0.25)
write.table(all.markers,
            file=paste0("./",sam.name,"/",sam.name,"_resolution0.5_total_marker_genes_tsne_",max(dim.use),"PC.txt"),
            sep="\t",quote = F,row.names = F)
save(experiment.aggregate,file = 'experiment.aggregate_before_anno.RData')
####手动细胞亚群注释####
library(dplyr)
library(openxlsx)
cellmarker.file <- "./multi/multi_resolution0.5_total_marker_genes_tsne_30PC.txt"
sampleid <- "multi"
#读入单细胞分析中输出的cell marker 基因文件
marker.gene <- read.table(cellmarker.file,header = T,stringsAsFactors = F,sep = "\t")
marker.gene.sig <- marker.gene %>% filter(as.numeric(p_val_adj) <= 0.05)

cat("Totally marker genes:",length(unique(marker.gene.sig$gene)))
table(marker.gene.sig$cluster)

#### 2. 定位数据库中存在的基因 ####
marker.gene.sig %>% filter(gene %in% cell.markers.tb$geneSymbol) -> marker.gene.sel
marker.gene.sel$cellMarker <- apply(marker.gene.sel,1,
                                    function(x)paste(unique(cell.markers.tb[cell.markers.tb$geneSymbol == x["gene"],1]),
                                                     collapse = ","))
cat("Totally find",length(unique(marker.gene.sel$gene)),"/",length(unique(marker.gene.sig$gene)),"genes in cellMarker db")

#### 3. 将marker基因与细胞类型结果写出到文件 ####
write.table(marker.gene.sel,file =paste0("./",sampleid,"/",sampleid,"_DEG_marker_cells_tsne_resolution0.5.txt"),row.names = T,col.names = T,sep = "\t",quote = F)

#样本的分组
meta1<-data.frame(matrix(nrow=length(experiment.aggregate@meta.data$orig.ident), ncol=2)) 
colnames(meta1)=c('Sample','Group')
meta1$Sample=experiment.aggregate@meta.data$orig.ident
unique(meta1$Sample)
### Group1 Tumor 为原发性肿瘤；Normal：正常
meta1[grep("NM01",meta1$Sample),]$Group="Healthy"
meta1[grep("NM01",meta1$Sample),]$Group="Healthy"
meta1[grep("NM01",meta1$Sample),]$Group="Healthy"
meta1[grep("DKD01",meta1$Sample),]$Group="DKD"
meta1[grep("DKD01",meta1$Sample),]$Group="DKD"
meta1[grep("DKD01",meta1$Sample),]$Group="DKD"
experiment.aggregate <- AddMetaData(experiment.aggregate, meta1$Sample,col.name = "Sample")
experiment.aggregate <- AddMetaData(experiment.aggregate, meta1$Group,col.name = "Group")

####根据基因表达手动匹配####
#Plasmacytoid dendritic cells
{
p3 = DotPlot(experiment.aggregate, features =  c(
  "IRF7"#,##T_Nk_cell
  #"CPA3","HPGDS","VWA5A",##Mast cells
  # "AIF1","MS4A6A","C1QA"##Kupffer cells
),cols = c( '#0f7ab0', '#ca443d', '#a51a49')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3
}
FeaturePlot(object = experiment.aggregate, 
            features=c("GATA3"),
        #reduction = "tsne", 
        #split.by  = "Group",
        label = T, 
        pt.size = 0.5) 
p3 = DotPlot(experiment.aggregate, features =  c(
  "CD1C","LYZ"
),cols = c( '#0f7ab0', '#ca443d', '#a51a49')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3

VlnPlot(experiment.aggregate, 
        features = c("CD1C","LYZ"),
        pt.size = 0,
        ncol = 2)

###通过标记基因及文献，可以人工确定各分类群的细胞类型，则可以如下手动添加细胞群名称###

###或者添加celltype
celltype=data.frame(ClusterID=0:31, ##0：15，多少个cluster，需要修改
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0,2,9,17,21),2]='Proximal tubule cell'
celltype[celltype$ClusterID %in% c(1,3,4),2]='Natural killer cell'
celltype[celltype$ClusterID %in% c(5,10,12),2]='Endothelial cell'
celltype[celltype$ClusterID %in% c(6,22),2]='Fibroblast'
celltype[celltype$ClusterID %in% c(7,28),2]='Mural cell'
celltype[celltype$ClusterID %in% c(8,13,16,23),2]='Distal convoluted tubule cell'
celltype[celltype$ClusterID %in% c(11,14,18,31,27),2]='Macrophage'
celltype[celltype$ClusterID %in% c(15,19,29,30),2]='B cell'
celltype[celltype$ClusterID %in% c(20,24),2]='Collecting duct cell'
celltype[celltype$ClusterID %in% c(25),2]='Mast cell'
celltype[celltype$ClusterID %in% c(26),2]='Granulocyte'

celltype 
table(celltype$celltype)
#先加一列celltype所有值为空
experiment.aggregate@meta.data$celltype = "NA"
###注释
for(i in 1:nrow(celltype)){
  experiment.aggregate@meta.data[which(experiment.aggregate@meta.data$seurat_clusters == celltype$ClusterID[i]),
                                 'celltype'] <- celltype$celltype[i]}
table(experiment.aggregate@meta.data$celltype)

#根据注释结果重命名每个细胞亚群
#Idents(experiment.aggregate) <- experiment.aggregate$seurat_clusters
new.cluster.ids <- c("Proximal tubule cell", "Natural killer cell", "Proximal tubule cell",#2
                     "Natural killer cell", "Natural killer cell", "Endothelial cell",#5
                     "Fibroblast", "Mural cell", "Distal convoluted tubule cell", #8
                     "Proximal tubule cell", "Endothelial cell", "Macrophage",#11
                     "Endothelial cell", "Distal convoluted tubule cell", "Macrophage", #14
                     "B cell", "Distal convoluted tubule cell", "Proximal tubule cell", #17
                     "Macrophage", "B cell", "Collecting duct cell", #20
                     "Proximal tubule cell", "Fibroblast", "Distal convoluted tubule cell", #23
                     "Collecting duct cel", "Mast cell", "Granulocyte", #26
                     "Macrophage", "Mural cell", "B cell", #29
                     "B cell", "Macrophage"#, #"Fibroblast" #32
                     )
names(new.cluster.ids) <- levels(experiment.aggregate)
experiment.aggregate <- RenameIdents(experiment.aggregate, new.cluster.ids)
Idents(experiment.aggregate)

save(experiment.aggregate,file = 'experiment.aggregate_afteranno_new.RData')

DimPlot(experiment.aggregate, reduction = "tsne", label = T,label.size = 5,cols = cluster_cols,group.by = "celltype", pt.size = 0.5) + NoLegend()
DimPlot(experiment.aggregate, reduction = "tsne", label = T,label.size = 56,group.by = "celltype",
        cols = cluster_cols,split.by = "Group", pt.size = 0.5) + NoLegend()
DimPlot(experiment.aggregate, reduction = "tsne", label = T,label.size = 6,group.by = "seurat_clusters",
        cols = cluster_cols,split.by = "orig.ident",ncol = 5, pt.size = 0.5) + NoLegend()


#按照数据来源分组展示细胞异同--画在多张图中
pdf(paste0("./",sam.name,"/afteranno_CellCluster-TSNEPlot_SamGroup_slipt_",max(dim.use),"PC.pdf"),width = 17,height = 5)
DimPlot(object = experiment.aggregate, 
        split.by ="Group",group.by = "celltype", 
        pt.size=0.5,reduction = "tsne",ncol = 2)
dev.off()

#可视化命名后的UMAP
pdf(file="after_anno_tsne_cell_type.pdf",width = 7,height = 7)
DimPlot(object = experiment.aggregate, reduction = "tsne", 
        group.by = "celltype",
        label = T, 
        pt.size = 0.5) 
#+ NoLegend()
dev.off()

#命名后的heatmap
pdf('after_anno_markers_heatmap_celltype.pdf', width = 16, height = 15)
top5 <- all.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC) 
DoHeatmap(experiment.aggregate, features = top5$gene, size = 3, angle = -50,hjust=0.8) + NoLegend()
dev.off()

####展示细胞群的差异基因####
# find markers
scRNA.markers <- FindAllMarkers(experiment.aggregate, only.pos = FALSE,
                                min.pct = 0.1,
                                logfc.threshold = 0)
write.table(scRNA.markers,file =paste0("cell_afternooo-DEG_genes-0.5.txt"),row.names = T,col.names = T,sep = "\t",quote = F)

p3 = DotPlot(experiment.aggregate, features =  c("RBMX"
),cols = c( '#0f7ab0', '#ca443d', '#a51a49')) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3

p3 = DotPlot(experiment.aggregate, group.by = "celltype",
             features =  c("FCER1G","FCER1A","FCGR1A","FCGR2A","FCGR3A","PLA2G4D","IL1R2"),
             cols = c( '#0f7ab0', '#ca443d', '#a51a49')) + 
  coord_flip() + #翻转
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p3

##选择基因
##小提琴图
pdf('genes_violin_plot_celltype.pdf', width = 10, height = 12)
VlnPlot(experiment.aggregate,
        features=c( "YTHDF2")#,ncol = 1
        )
dev.off()

FeaturePlot(experiment.aggregate, features = "RBM15")
####细胞群之间的差异分析####
diff_exp_results <- FindMarkers(experiment.aggregate, ident.1 = "Fibroblast", ident.2 = "T_cell")
write.csv(diff_exp_results,file = "fib_t_degs.csv")

####按细胞类型展示不同样本细胞比例####
# 将样本、细胞类型、细胞个数制成表格
sample_table4 <- as.data.frame(table(scRNA@meta.data$seurat_clusters,scRNA@meta.data$Group))
# 定义列名
names(sample_table4 ) <- c("seurat_clusters","Group","CellNumber")
# 指定颜色
sample_colors <- c("#DC050C",  "#1965B0","#FB8072", "#7BAFDE", "#882E72",
                   "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                   "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                   "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                   "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                   "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00","#D9A86C","#F2C450")
# 横轴为细胞类型,纵轴为样本比例
p6 <-ggplot(sample_table4,aes(x=seurat_clusters,weight=CellNumber,fill=Group))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
p6

# 将样本、细胞类型、细胞个数制成表格
sample_table2 <- as.data.frame(table(scRNA@meta.data$celltype,scRNA@meta.data$orig.ident))
# 定义列名
names(sample_table2 ) <- c("celltype","Samples","CellNumber")
# 指定颜色
sample_colors <- c("#DC050C", "#FB8072", "#1965B0", "#7BAFDE", "#882E72",
                   "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                   "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                   "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                   "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                   "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
# 横轴为细胞类型,纵轴为样本比例
p6 <-ggplot(sample_table2,aes(x=celltype,weight=CellNumber,fill=Samples))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
p6

#横轴为细胞类型，纵轴为不同样本来源细胞绝对数量：
p7 <-ggplot(sample_table2,aes(x=celltype,weight=CellNumber,fill=Samples))+
  geom_bar(position="dodge",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Cell number")+RotatedAxis()
p7 

##横轴为不同样本来源细胞比例，纵轴为细胞类型：
p8 <-ggplot(sample_table2,aes(x=celltype,weight=CellNumber,fill=Samples))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()+coord_flip() 
p8


##按细胞类型展示不同疾病分组比例
# 将样本、细胞类型、细胞个数制成表格
sample_table3 <- as.data.frame(table(experiment.aggregate@meta.data$celltype,
                                     experiment.aggregate@meta.data$Group))
# 定义列名
names(sample_table3 ) <- c("celltype","Group","CellNumber")
# 指定颜色
sample_colors_group <- c("#FB8072","#7BAFDE","#DC050C", "#1965B0", "#882E72",
                         "#B17BA6", "#FF7F00", "#FDB462", "#E7298A", "#E78AC3",
                         "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D",
                         "#E6AB02", "#7570B3", "#BEAED4", "#666666", "#999999",
                         "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000",
                         "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")
# 横轴为细胞类型,纵轴为样本比例
p6 <-ggplot(sample_table3,aes(x=celltype,weight=CellNumber,fill=Group))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors_group) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()
p6

#横轴为细胞类型，纵轴为不同样本来源细胞绝对数量：
p7 <-ggplot(sample_table3,aes(x=celltype,weight=CellNumber,fill=Group))+
  geom_bar(position="dodge",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors_group) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Cell number")+RotatedAxis()
p7 

##横轴为不同样本来源细胞比例，纵轴为细胞类型：
p8 <-ggplot(sample_table3,aes(x=celltype,weight=CellNumber,fill=Group))+
  geom_bar(position="fill",width = 0.7,size = 0.5,colour = '#222222')+
  scale_fill_manual(values = sample_colors_group) + 
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16)
  )+labs(y="Percentage")+RotatedAxis()+coord_flip() 
p8




####GSEA&GSVA####
#加载需要的R包
library(Seurat)
library(msigdbr)
library(GSVA)
library(tidyverse)
library(clusterProfiler)
library(patchwork)
library(fgsea)
library(dplyr)
library(ggplot2)
library(enrichplot)
setwd("E:/R_2024/1029_skin/GSEA_GSVA")
#加载单个样本单细胞数据分析结果
#load("scRNA.RData")
#挑选B细胞相对于CD8T细胞特意性高表达的marker基因
#B_vs_CD8T=FindMarkers(scRNA, ident.1="Mon",ident.2="KC",only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
#制作geneList
#值为log2FC
geneList= diff_exp_results$avg_log2FC 
#name为基因名字
names(geneList)= rownames(diff_exp_results)
#按FC从大到小排序
geneList=sort(geneList,decreasing = T)
head(geneList)
#从GSEA官网下载GSEA分析需要的基因集
#http://www.gsea-msigdb.org/gsea/index.jsp
#下载免疫相关的基因集，c7: immunologic signature gene sets

####进行GSEA_BP分析####
#读取gmt文件中的pathway信息
gmtfile ='GSEA_data/h.all.v2024.1.Hs.symbols.gmt'
pathway<-read.gmt(gmtfile)
y <- GSEA(geneList,TERM2GENE =pathway,pvalueCutoff = 0.05)
#保存GSEA分析结果
write.csv(file="fib_vs_t_GSEA_hallmark_result.csv",data.frame(y))
#气泡图展示显著富集的前20条通路
pdf(file="GSEA_dotplot_hallmark.pdf",width=10,he=15)
dotplot(y,showCategory=20)
dev.off()
#绘制具体通路的GSEA图
library(enrichplot)
pdf(file="fib_vs_t_hallmark1.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'WP_CCL18_SIGNALING',pvalue_table=T)
dev.off()
pdf(file="fib_vs_t_hallmark6.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'WP_GLYCOLYSIS_AND_GLUCONEOGENESIS',pvalue_table=T) 
dev.off()

####进行GSEA_CC分析####
gmtfile ='GSEA_data/m5.go.cc.v2024.1.Mm.symbols.gmt'
pathway<-read.gmt(gmtfile)
y <- GSEA(geneList,TERM2GENE =pathway,pvalueCutoff = 0.99)
#保存GSEA分析结果
write.csv(file="Mon_vs_KC_GSEA_cc_result.csv",data.frame(y))
#气泡图展示显著富集的前20条通路
pdf(file="GSEA_dotplot_cc.pdf",width=10,he=15)
dotplot(y,showCategory=20)
dev.off()
#绘制具体通路的GSEA图
library(enrichplot)
pdf(file="Mon_VS_KC_cc1.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'HALLMARK_APOPTOSIS',pvalue_table=T)
dev.off()
pdf(file="Mon_VS_KC_cc2.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'GOCC_LYTIC_VACUOLE',pvalue_table=T) 
dev.off()

####进行GSEA_CC分析####
gmtfile ='GSEA_data/m5.go.mf.v2024.1.Mm.symbols.gmt'
pathway<-read.gmt(gmtfile)
y <- GSEA(geneList,TERM2GENE =pathway,pvalueCutoff = 0.99)
#保存GSEA分析结果
write.csv(file="Mon_vs_KC_GSEA_mf_result.csv",data.frame(y))
#气泡图展示显著富集的前20条通路
pdf(file="GSEA_dotplot_mf.pdf",width=10,he=15)
dotplot(y,showCategory=20)
dev.off()
#绘制具体通路的GSEA图
library(enrichplot)
pdf(file="Mon_VS_KC_mf1.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'GOMF_CALCIUM_DEPENDENT_PROTEIN_BINDING',pvalue_table=T)
dev.off()
pdf(file="Mon_VS_KC_mf2.pdf",width=12,he=11)
gseaplot2(y,geneSetID = 'GOMF_SIGNALING_RECEPTOR_BINDING',pvalue_table=T) 
dev.off()

#GSVA
##########################################
#获取细胞类型
Idents(experiment.aggregate)
#计算每个基因在每个细胞亚群中的平均表达值
expr <- AverageExpression(experiment.aggregate, assays = "RNA", slot = "data")[[1]]
View(expr)
#选取非零基因
expr <- expr[rowSums(expr)>0.5,] 
#转换成矩阵
expr <- as.matrix(expr)

#从GSEA官网下载GSEA分析需要的基因集
#http://www.gsea-msigdb.org/gsea/index.jsp
setwd("E:/R_2024/1029_skin/GSEA_GSVA")
#下载免疫相关的基因集，h: hallmark gene sets
gmtfile ='c2.cp.kegg_legacy.v2024.1.Hs.symbols.gmt'

#读取gmt文件中的pathway信息
pathway<-read.gmt(gmtfile)[,c(2,1)]
#去堆叠,转换成list
genesets=unstack(pathway)

#进行GSVA分析
#gsva.res <- gsva(expr, genesets, method="ssgsea", useNames = TRUE) 
gsva.res <- gsva(expr, genesets, method="ssgsea") 
#保存GSVA分析结果
gsva.df <- data.frame(Genesets=rownames(gsva.res), gsva.res, check.names = F)
write.csv(gsva.df, "gsva_res_wikipathways.csv", row.names = F)

#绘制热图
library(pheatmap)
pdf(file="GSVA_heatmap_reactome.pdf",width=12,height=8)
pheatmap(gsva.res[151:188,], show_colnames = T, scale = "row")
dev.off()

####score####
#BiocManager::install("UCell")
library(UCell)
#做基因集评分
genes <- list(c("METTL3", "METTL14", "RBM15",
                "RBM15B", "WTAP", "KIAA1429",
                "CBLL1", "ZC3H13","ALKBH5", 
                "FTO","YTHDC1","YTHDC2",
                "YTHDF1", "YTHDF2", "YTHDF3", 
                "IGF2BP1","HNRNPA2B1", "HNRNPC", 
                "FMR1", "LRPPRC", "ELAVL1"))
names(genes) <- 'gene'
gene_score <- AddModuleScore_UCell(experiment.aggregate,features=genes,name="_score")
#提取数据
library(ggpubr)
df<- FetchData(gene_score,vars = c("orig.ident","gene_score"))
df$orig.ident <- factor(df$orig.ident,levels = c("HC","EEC","AEH"))#设置顺序
#设置比较组
my_comparisons1 <- list(c("HC", "EEC"))
my_comparisons2 <- list(c("AEH", "EEC"))

####细胞通讯####
save(experiment.aggregate,file = 'experiment.aggregate_before_cellchat_new.RData')
setwd("E:/R_2024/1029_skin/cellchat")
library(CellChat)
library(tidyverse)
library(Seurat)
options(stringsAsFactors = FALSE)
cellchat <- createCellChat(object = experiment.aggregate,                           
                           meta = experiment.aggregate@meta.data,                           
                           group.by = "celltype")
cellchat

CellChatDB <- CellChatDB.human  
showDatabaseCategory(CellChatDB)
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat@DB <- CellChatDB
#注1：如果想用全部的用于cellchat分析，不进行subsetDB，直接指定cellchat@DB <- CellChatDB 即可。

cellchat <- subsetData(cellchat) 
future::plan("multiprocess", workers =10) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat@data.project[1:4,1:4]

###二 推断cell-cell communication network
##1.推断细胞通讯网络
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
cellchat <- filterCommunication(cellchat, min.cells = 10)
##2.提取 保存结果
#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "cell-cell_communications.all.csv")
#access the the inferred communications at the level of signaling pathways
df.net1 <- subsetCommunication(cellchat,slot.name = "netP")
#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
levels(cellchat@idents)
df.net2 <- subsetCommunication(cellchat, sources.use = c("Fibroblast"), targets.use = c("T_cell")) 
#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
##3.计算cell-cell communication
#计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- computeCommunProbPathway(cellchat)
#计算整合的细胞类型之间通信结果
cellchat <- aggregateNet(cellchat)

###CellChat 可视化
##1.celltype之间通讯结果
#1）根据使用netVisual_circle显示任意两个celltype之间的通讯次数（左）或总通讯强度(右)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, 
                 weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize,
                 weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

#2)分别展示根据celltype的个数灵活调整mfrow = c(2,3) 参数。像绘制count维度的，只需要修改下mat <- cellchat@net$count即可。
mat <- cellchat@net$weight
par(mfrow = c(2,3), xpd=TRUE)for (i in 1:nrow(mat)) {  
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))  
  mat2[i, ] <- mat[i, ]  
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
  }
##2.单个信号通路可视化
#根据cellchat@netP$pathways展示当前有哪些通路结果
#1层级图绘制层级图的话 ,需要指定layout为hierarchy ，当前版本默认下出来的是circle图。
cellchat@netP$pathways
pathways.show <- c("VEGF")  
levels(cellchat@idents)  
vertex.receiver = c(1,2,4,5) 
netVisual_aggregate(cellchat, signaling = pathways.show,                      
                    vertex.receiver = vertex.receiver,layout = "hierarchy")
#左图中间的Target是vertex.receiver选定的细胞类型，右图是除vertex.receiver选中之外的另外的细胞类型
#2）和弦图和热图
#Chord diagram
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
#Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
##3.绘制配体受体气泡图
#1）指定受体-配体细胞类型
#绘制指定受体-配体细胞类型中的全部配体受体结果的气泡图,通过sources.use 和 targets.use指定。
#指定受体-配体细胞类型
netVisual_bubble(cellchat, sources.use = c(4),
 targets.use = c(1,2,3,5,6,7,8,9), remove.isolate = FALSE)
#2）指定受体-配体细胞类型 且 指定通路
#指定TGFb和SPP1两个信号通路
netVisual_bubble(cellchat, sources.use = c(4), targets.use = c(9),                  
signaling = c("CXCL","MIF","CCL","IL1","TNF","IL6","VEGF"), remove.isolate = FALSE)
#3）某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
p1 = plotGeneExpression(cellchat, signaling = "VEGF")
p1
#p2 = plotGeneExpression(cellchat, signaling = "CXCL", type = "dot")


#####细胞通讯2####
#导入需要的R包
setwd("D:\0502\0526_DKD\GSE209781\cellchat")
####cellchat2#####
#导入需要的R包
library(CellChat)
library(patchwork)
#devtools::install_github("tidyverse/magrittr")
library(magrittr)
library(tidyverse)
library(Seurat)

#读取示例数据
#该数据为前期保存的单细胞数据分析结果
#setwd("D:/Resoure/DR/mouse/0718/1019/cellchat2")
immune.combined<- experiment.aggregate

#Chellchat多样本分析
## STIM组样本chellchat分析
stim.object <- subset(immune.combined,Group1=="DKD")
stim.data.input <- GetAssayData(stim.object, assay = "RNA", slot = "data")
stim.meta <- stim.object@meta.data[,c("celltype", "Group1")]
stim.meta$CellType %<>% as.vector(.)

#创建CellChat对象
stim.cellchat <- createCellChat(object = stim.data.input)
stim.cellchat <- addMeta(stim.cellchat, meta = stim.meta)
stim.cellchat <- setIdent(stim.cellchat, ident.use = "celltype")
cellchat <- createCellChat(object = stim.data.input, meta = stim.meta, group.by = "Group1")

# CellChat提供的人的配受体数据库
stim.cellchat@DB <- CellChatDB.human
#stim.cellchat@DB <- CellChatDB.mouse
stim.cellchat <- subsetData(stim.cellchat)
#使用多线程进行计算

future::plan("multisession", workers = 10)

#用于细胞-细胞通信分析的表达数据预处理
options(future.globals.maxSize = 8000 * 1024^2)
stim.cellchat <- identifyOverExpressedGenes(stim.cellchat)
stim.cellchat <- identifyOverExpressedInteractions(stim.cellchat)
stim.cellchat <- projectData(stim.cellchat, PPI.human)

#在信号通路水平推断细胞间通信
stim.cellchat <- computeCommunProb(stim.cellchat)

#计算聚合的细胞间的通信网络
stim.cellchat <- aggregateNet(stim.cellchat)

#保存cellchat分析结果
saveRDS(stim.cellchat,"stim.cellchat.rds")

## CTRL组样本Cellchat分析
ctrl.object <- subset(immune.combined,Group1=="Healthy")
ctrl.data.input <- GetAssayData(ctrl.object, assay = "RNA", slot = "data")
ctrl.meta = ctrl.object@meta.data[,c("celltype", "Group1")]
ctrl.meta$CellType %<>% as.vector(.)
ctrl.cellchat <- createCellChat(object = ctrl.data.input)
ctrl.cellchat <- addMeta(ctrl.cellchat, meta = ctrl.meta)
ctrl.cellchat <- setIdent(ctrl.cellchat, ident.use = "celltype")
ctrl.cellchat@DB <- CellChatDB.human
ctrl.cellchat <- subsetData(ctrl.cellchat)
future::plan("multisession", workers = 10)

#用于细胞-细胞通信分析的表达数据预处理
ctrl.cellchat <- identifyOverExpressedGenes(ctrl.cellchat)
ctrl.cellchat <- identifyOverExpressedInteractions(ctrl.cellchat)
ctrl.cellchat <- projectData(ctrl.cellchat, PPI.human)
ctrl.cellchat <- computeCommunProb(ctrl.cellchat)
ctrl.cellchat <- computeCommunProbPathway(ctrl.cellchat)

#计算聚合的细胞间的通信网络
ctrl.cellchat <- aggregateNet(ctrl.cellchat)

#保存cellchat分析结果
saveRDS(ctrl.cellchat,"ctrl.cellchat.rds")

#加载两个数据集的Cellchat对象，然后合并在一起
object.list <- list(Healthy = ctrl.cellchat, DKD = stim.cellchat)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

#进行互作次数和互作强度分析
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
pdf("compareInteractions.pdf",height = 5,width = 8)

gg1 + gg2
dev.off()
#注：两个分组的互作次数和互作强度柱状图
#不同细胞群相互作用次数和作用强度差异
#网络图绘制
pdf("netVisual_diffInteraction.pdf",height = 10,width = 10)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
dev.off()

##热图绘制
gg1 <- netVisual_heatmap(cellchat)
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
pdf("netVisual_heatmap.pdf",height = 5,width = 10)
gg1 + gg2
dev.off()

#不同数据集的互作网络图

pdf("netVisual_diffInteraction.pdf",height = 10,width = 10)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, 
                   label.edge= F, edge.weight.max = weight.max[2], 
                   edge.width.max = 12, 
                   title.name = paste0("Number of interactions – ", names(object.list)[i]))
  
}
dev.off()

levels(cellchat@idents)  
#识别上调和下调的受体-配体对
pdf("netVisual_bubble.pdf",height = 10,width = 10)
netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10), comparison = c(1, 2), angle.x = 45)
dev.off()

#气泡图展示所有配体受体对的差异
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10),
                      comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble.pdf", p, width = 12, height = 8)

#气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use = 11, targets.use = c(1:10), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in DKD", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat,sources.use = 11, targets.use = c(1:10), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in DKD", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated.pdf", pc, width = 15, height = 7.5)

#气泡图展示所有配体受体对的差异
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use =  c(1:10), targets.use = 11,
                      comparison = c(1, 2), angle.x = 45)
ggsave("Compare_LR_bubble1.pdf", p, width = 16, height = 8)

#气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use =  c(1:10), targets.use = 11, comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in DKD", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat,sources.use =  c(1:10), targets.use = 11, comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in DKD", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Compare_LR_regulated1.pdf", pc, width = 17, height = 7.5)

