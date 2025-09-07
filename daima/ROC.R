
library(pROC)               
expFile="GSE30122.txt"    
geneFile="interGenes.txt"     
setwd("D:/0502/0526_DKD/GSE209781/RNA_seq/07.rna-seq/test/GSE30122") 
rt=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
y=gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y=ifelse(y=="Con", 0, 1)

geneRT=read.table(geneFile, header=F, sep="\t", check.names=F)

for(x in as.vector(geneRT[,1])){
	roc1=roc(y, as.numeric(rt[x,]))
	ci1=ci.auc(roc1, method="bootstrap")
	ciVec=as.numeric(ci1)
	pdf(file=paste0("ROC.",x,".pdf"), width=5, height=5)
	plot(roc1, print.auc=TRUE, col="red", legacy.axes=T, main=x)
	text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
	dev.off()
}

