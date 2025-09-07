library(limma)
library(pheatmap)
library(ggplot2)

inputFile="geoMatrix.txt"     #表达数据文件
conFile="s1.txt"               #对照组的样品文件
treatFile="s2.txt"             #实验组的样品文件
logFCfilter=1                  #logFC过滤条件(logFC=0.585,差异倍数1.5倍;logFC=1,差异2倍;logFC=2,差异4倍)
adj.P.Val.Filter=0.05          #矫正后p值的过滤条件

geoID="gse37171"               #GEO数据库研究的id
conName="con"              #对照组名称
treatName="treat"          #实验组名称
setwd("F:\\data\\人\\result_06\\RNA_seq\\rna-seq\\01.data")      #设置工作目录

#读取输入文件，并对输入文件整理
rt=read.table(inputFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
rt=data[rowMeans(data)>0,]

#如果数据没有取log2, 会对数据自动取log2
qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}
data=normalizeBetweenArrays(rt)

#读取样品信息的文件
sample1=read.table(conFile, header=F, sep="\t", check.names=F)
sample2=read.table(treatFile, header=F, sep="\t", check.names=F)
#sampleName1=gsub("^ | $", "", as.vector(sample1[,1]))
sampleName1= as.vector(sample1[,1])
sampleName2=as.vector(sample2[,1])
#sampleName2=gsub("^ | $", "", as.vector(sample2[,1]))
conData=data[,sampleName1]
treatData=data[,sampleName2]
data=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#差异分析
Type=c(rep("con",conNum), rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(data,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
allDiffOut=rbind(id=colnames(allDiff),allDiff)
write.table(allDiffOut,file=paste0(geoID,".all.txt"),sep="\t",quote=F,col.names=F)

#输出矫正后的表达量
Type=c(rep(conName,conNum),rep(treatName,treatNum))
outData=rbind(id=paste0(colnames(data),"_",Type),data)
write.table(outData,file=paste0(geoID,".normalize.txt"),sep="\t",quote=F,col.names=F)

#输出显著的差异基因
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & P.Value < adj.P.Val.Filter )), ]
diffSigOut=rbind(id=colnames(diffSig),diffSig)
write.table(diffSigOut,file=paste0(geoID,".diff.txt"),sep="\t",quote=F,col.names=F)

#绘制差异基因热图
geneNum=50    #设置基因的数目
diffSig=diffSig[order(as.numeric(as.vector(diffSig$logFC))),]
diffGeneName=as.vector(rownames(diffSig))
diffLength=length(diffGeneName)
hmGene=c()
if(diffLength>(2*geneNum)){
  hmGene=diffGeneName[c(1:geneNum,(diffLength-geneNum+1):diffLength)]
}else{
  hmGene=diffGeneName
}
hmExp=data[hmGene,]
#定义注释文件
Type=c(rep(conName,conNum),rep(treatName,treatNum))
Type=factor(Type, levels=c(conName, treatName))
names(Type)=colnames(data)
Type=as.data.frame(Type)
ann_colors = list(Type = c(con="#2a83a2", treat="#bb5548"))

windows()
#绘制热图
pdf(file=paste0(geoID,".heatmap.pdf"), width=9, height=7)
pheatmap(hmExp, 
         annotation=Type, 
         annotation_colors=ann_colors,
         color = colorRampPalette(c(rep("#2a83a2",15), "white",rep("#bb5548",15)))(100),
         cluster_cols =F,
         show_colnames = F,
         scale="row",
         gaps_col=as.vector(cumsum(table(Type))),
         fontsize = 10,
         fontsize_row=6,
         fontfamily="serif",
         #fontface="italic",
         #fontface="bold",
         #font=2,
         #fontfamily= "Times New Roman",
         fontsize_col=7)
dev.off()

#定义显著性
allDiff$logFC[allDiff$logFC>20]=20
allDiff$logFC[allDiff$logFC< -20]=-20
Significant=ifelse((allDiff$P.Value<adj.P.Val.Filter & abs(allDiff$logFC)>logFCfilter), ifelse(allDiff$logFC>logFCfilter,"Up","Down"), "Not")
#绘制火山图
p = ggplot(allDiff, aes(logFC, -log10(P.Value))) +
  geom_point(aes(col = Significant)) +
  scale_color_manual(values = c("#2a83a2", "grey", "#bb5548")) +
  labs(title = " ") +
  theme(plot.title = element_text(size = 10, hjust = 0.5, family = "Times New Roman", face = "bold")) +
  theme_bw() +
  theme(text = element_text(size = 10))  # 设置所有文本字体大小为10pt

# 保存为图片
 pdf(paste0(geoID, ".vol.pdf"), width = 10, height = 5)  # 宽度设为10，高度保持5
print(p)
dev.off()
pdf("F:/data/人/result_06/RNA_seq/rna-seq/01.data/volcano_plot.pdf", width = 10, height = 5)
print(p)
dev.off()