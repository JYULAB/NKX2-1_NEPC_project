setwd("D:\\bigdata\\dimelo\\2day_vs_4w\\delta_02_span200\\PDX_cluster")

library(methylKit)
library(genomation)

file.list=list("D2_FOXA2_5mc_DMR_region_only.txt", "D14_FOXA2_5mc_DMR_region_only.txt", "D28_FOXA2_5mc_DMR_region_only.txt", "CRPC_35CR_5mc_DMR_region_only.txt", "CRPC_70CR_5mc_DMR_region_only.txt", "CRPC_77CR_5mc_DMR_region_only.txt", "CRPC_147_5mc_DMR_region_only.txt", "NEPC_93_5mc_DMR_region_only.txt", "NEPC_145_1_5mc_DMR_region_only.txt", "NEPC_145_2_5mc_DMR_region_only.txt", "NEPC_173_1_5mc_DMR_region_only.txt", "NCI_H660_5mc_DMR_region_only.txt")
myobj=methRead(file.list, sample.id=list("D2_FOXA2", "D14_FOXA2", "D28_FOXA2", "CRPC_35CR", "CRPC_70CR", "CRPC_77CR", "CRPC_147", "NEPC_93", "NEPC_145.1", "NEPC_145.2", "NEPC_173.1", "NCI_H660"), assembly="hg38", treatment=c(0,0,1,0,0,0,0,1,1,1,1,1), context="CpG", pipeline=list(fraction=FALSE,chr.col=1,start.col=2,end.col=3,coverage.col=10,strand.col=6,freqC.col=11), header=FALSE, mincov = 10)

meth=unite(myobj, destrand=TRUE)
meth.min=unite(myobj,min.per.group=1L)
clusterSamples(meth, dist="correlation", method="average", sd.threshold=0.5,filterByQuantile=TRUE,plot=TRUE)

PCASamples(meth,screeplot=FALSE, adj.lim=c(0.0004,0.1),
scale=TRUE,center=TRUE,comp=c(1,2),transpose=FALSE,sd.filter=TRUE,
sd.threshold=0.5,filterByQuantile=TRUE,obj.return=FALSE)

getCorrelation(meth,plot=TRUE)






setwd("D:\\bigdata\\dimelo\\2day_vs_4w\\delta_02_span200")
library(ggplot2)

data=read.table("pca_RNA_seq_felix_foxa2_pdx_input.txt",header=T,sep="\t",row.names=1)  
data=t(as.matrix(data))   
data.class <- rownames(data)
data.pca <- prcomp(data, scale. = TRUE)   
pca.sum=summary(data.pca)


#pca 2d plot
#PDX_173_1 is far away from other ARneg_NEpos, remove it
#group=c(rep("Felix_NEPC",5),rep("Felix_CRPC",93),rep("FOXA2_D2",3),rep("FOXA2_D28",3),rep("PDX_49",2),rep("PDX_93",2),rep("PDX_145_1",2),rep("PDX_145_2",2),rep("PDX_93",2),rep("PDX_23",3),rep("PDX_35CR",3),rep("PDX_70CR",2),rep("PDX_73CR",2),rep("PDX_78CR",2),rep("PDX_81CR",2),rep("PDX_96CR",3),rep("PDX_105CR",2),rep("PDX_136CR",2),rep("PDX_147CR",2)) 
group=c(rep("Felix_NEPC",5),rep("Felix_CRPC",93),rep("FOXA2_D2",3),rep("FOXA2_D7",3),rep("FOXA2_D14",3),rep("FOXA2_D21",3),rep("FOXA2_D28",3),rep("PDX_ARneg_NEpos",8),rep("PDX_ARpos_NEneg",23)) 
pcaPredict=predict(data.pca)
PCA = data.frame(PCA1 = pcaPredict[,1], PCA2 = pcaPredict[,2],group=group)
PCA.mean=aggregate(PCA[,1:2],list(group=PCA$group),mean)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 2, npoints = 600) 
  {
    theta <- (0:npoints) * 2 * pi/npoints
    Circle <- cbind(cos(theta), sin(theta))
    t(center + scale * t(Circle %*% chol(cov)))
  }
df_ell <- data.frame()
for(g in levels(PCA$group)){
  df_ell <- rbind(df_ell, cbind(as.data.frame(with(PCA[PCA$group==g,],
                  veganCovEllipse(cov.wt(cbind(PCA1,PCA2),
                  wt=rep(1/length(PCA1),length(PCA1)))$cov,
                  center=c(mean(PCA1),mean(PCA2))))),group=g))
}

pdf(file="pca_RNA_seq_felix_foxa2_pdx.pdf")
ggplot(data = PCA, aes(PCA1, PCA2)) + geom_point(aes(color = group)) +
    geom_path(data=df_ell, aes(x=PCA1, y=PCA2,colour=group), size=0.6, linetype=2)+
    annotate("text",x=PCA.mean$PCA1,y=PCA.mean$PCA2,label=PCA.mean$group)+
    theme_bw()+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()

