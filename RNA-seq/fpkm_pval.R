# Making plot of average FPKM and p-values

#Load libraries
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(ggprism)

#Change working directory
setwd("E:/Yu Lab/rnaseq/scatter")

#Read in data
all <- read.delim("CRPC and NEPC patient 2019 JCI.txt", header=TRUE, sep="\t", na.strings="NA")

#Isolate rows of FOXA2 and other gene
df <- all %>% filter(grepl('FOXA2|FOXA1', X))
df <- df[-3,]
row.names(df) <- df[,1]
df <- df[,-1]
df <- as.data.frame(t(df))
df <- data.frame(sapply(df, as.numeric), row.names=rownames(df))

#Identify AR and NE character of each sample
metadata <- read.delim("meta_pt.txt", header=TRUE, sep="\t", na.strings="NA")
row.names(metadata) <- metadata[,1]
fpkm <- merge(df, metadata, by=0)

#fpkm <- fpkm[,-c(1,4)]
rownames(fpkm) <- fpkm[,1]
fpkm <- fpkm[,-1]
fpkm <- fpkm[,-3]

fpkm[,5]<- (fpkm[,1] - mean(fpkm[,1]))/sd(fpkm[,1])
fpkm[,6]<- (fpkm[,2] - mean(fpkm[,2]))/sd(fpkm[,2])

fpkm[,7]<- fpkm[,1]+1
fpkm[,8]<- fpkm[,2]+1
colnames(fpkm)<- c("rlog_FOXA2","rlog_FOXA1","cond", "PT","zscore_FOXA2","zscore_FOXA1","rlog_FOXA2+1","rlog_NKX2.1+1")


#write.table(fpkm, file = "rlog_zscores.txt", append = FALSE, quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE, na = "")


#### Wilcoxon Test for p-values ####

f.stat.test <- fpkm %>%
  wilcox_test(zscore_FOXA2 ~ cond, p.adjust.method = "none")
f.stat.test <- f.stat.test %>% filter(p.adj < 0.05)
f_pval <- as.data.frame(f.stat.test[c(2,3,8,9)])
#f_pval$xmin <- c(30,30,70,70,75)
#f_pval$xmax <- c(85,95,85,95,95)
f_pval$xmin <- c(5,5,18,27)
f_pval$xmax <- c(18,27,11,11)
f_pval$p.adj <- signif(f_pval$p.adj, digits=3)
f_pval <- f_pval[order(abs(f_pval$xmax-f_pval$xmin)),]
f_pval$ypos <- c(1,1.5,2,2.5)
write.table(f_pval, "labrecque_foxa2_pval.txt", sep=" \t")

n.stat.test <- fpkm %>%
  wilcox_test(zscore_FOXA1 ~ cond, p.adjust.method = "none")
n.stat.test <- n.stat.test %>% filter(p.adj < 0.05)
n_pval <- as.data.frame(n.stat.test[c(2,3,8,9)])
#n_pval$xmin <- c(30,30,70,70,75,85)
#n_pval$xmax <- c(75,90,75,90,90,90)
n_pval$xmin <- c(5,5)
n_pval$xmax <- c(27,11)
n_pval$p.adj <- signif(n_pval$p.adj, digits=3)
n_pval <- n_pval[order(abs(n_pval$xmax-n_pval$xmin)),]
n_pval$ypos <- c(11,12)
write.table(n_pval, "labrecque_foxa1_pval.txt", sep=" \t")


####Plot####
# need to make conditions to get the correct order of plot
#fpkm$cond <- factor(fpkm$cond, levels=c("AR+/NE+", "AR+/NE-", "AR-/NE+", "ARlow/NE-","AR-/NE-"))

#fpkm <- fpkm[order(fpkm$cond, fpkm$log2_FOXA2), ]

#Correlation Plot
require('RColorBrewer') #color palette for the scatter plot
pick.col <- c("#FF66CC","#FF6600", "#6699FF","#9933FF")
"#000000""#CC9900""#006633"
colII <- colorRampPalette(pick.col)(4)[factor(4, levels = metadata)]

ggplot(fpkm, aes(x=rlog_FOXA2, y=rlog_NKX2.1)) + 
  geom_point(aes(color = cond),size=2) +
  theme_minimal() +
  scale_color_manual(values = pick.col) + 
  labs(x="FOXA2 rlog(CPM)", y="NKX2.1 rlog(CPM)")+
  geom_smooth(method='lm', formula= y~x, color="black")+
  stat_cor(method = "pearson", label.x = 3, label.y = 3)
#geom_text(aes(label=select), size=4, vjust=1.5, hjust=0.5)



#fpkm$cond <- factor(fpkm$cond, levels = c("AR+/NE-", "AR+/NE+", "ARlow/NE-", "AR-/NE-", "AR-/NE+"))
fpkm$cond <- factor(fpkm$cond, levels = c("CRPC-AR","CRPC-SCL", "CRPC-WNT", "CRPC-NE"))

#FOXA2
pick.col <- c("#F8766D","#FF6600", "#00B0F6","#FB61D7")
fpkm$PT <- factor(fpkm$PT, levels = fpkm$PT[order(fpkm$rlog_FOXA2)])
cond <- fpkm$cond
n <- fpkm %>% group_by(cond) %>% tally()
n$grouping <- sprintf("%s, n= %s", n$cond, n$n)
# n$colors <- c("coral", "darkolivegreen4","seagreen3","turqouise3","darkorchid3")
  
f <- ggplot(fpkm, 
            aes(x=interaction(PT, cond),
                y=zscore_FOXA2)) +
  geom_col(aes(fill=cond), position="dodge")+
  xlab("condition")+
  ylab("zscore(rlog(CPM))")+
  #ylim(-0.5, 7.25)+
  scale_fill_manual(values = pick.col) +
  theme( panel.grid.major.x = element_blank(),
         axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

f

f + geom_bracket(
  xmin = f_pval$xmin, xmax = f_pval$xmax,
  y.position = f_pval$ypos, 
  label = format(f_pval$p.adj.signif, scientific = TRUE),
  tip.length = 0.01,
  vjust = 0.5,
  hjust = 0.3,
  label.size =6
) 



#NKX2-1
fpkm$PT <- factor(fpkm$PT, levels = fpkm$PT[order(fpkm$rlog_NKX2.1)])

nk <- ggplot(fpkm, 
            aes(x=interaction(PT, cond),
                y=rlog_NKX2.1)) +
  geom_col(aes(fill=cond), position="dodge")+
  xlab("condition")+
  ylab("rlog(CPM)")+
  ylim(0,15)+
  scale_fill_discrete(labels = n$grouping)+
  theme( panel.grid.major.x = element_blank(),
         axis.ticks.x = element_blank(),
         axis.text.x = element_blank())

nk

nk + geom_bracket(
  xmin = n_pval$xmin, xmax = n_pval$xmax,
  y.position = n_pval$ypos, 
  label = format(n_pval$p.adj.signif, scientific = TRUE),
  tip.length = 0.01,
  vjust = 0.5,
  hjust = 0.3,
  label.size =6
)

table_p <- f_pval[,1:4]

write.table(table_p,"E://Yu Lab//rnaseq//scatter//log2FPKM_FOXA2_pvals.txt", append = FALSE, quote = TRUE, sep = " ")

