library(RSQLite)
library(dplyr)
library("ggpubr")
library("tidyr")

setwd("D:\\bigdata\\dimelo\\new_code")    

gene_name = 'NKX2-1'

DB <- dbConnect(SQLite(), dbname = "D:\\pdb\\database\\METHYLATION.db")
NEPC_sample <- c('SRR3146982','SRR3146983','SRR3146984','SRR3146993','SRR3146995','SRR3146996','SRR3146997','SRR3146998','SRR3147004','SRR3147007')
CRPC_sample <- c('SRR3146985','SRR3146986','SRR3146987','SRR3146988','SRR3146989','SRR3146990','SRR3146991','SRR3146992','SRR3146994','SRR3146999','SRR3147000','SRR3147001','SRR3147002','SRR3147003','SRR3147005','SRR3147006','SRR3147008','SRR3147009')
statement2 <- paste0("SELECT n.id, n.name, disToIsland, readStart pos, percent frac, sampleName from CPG_ISLAND_BELTRAN_RRBS m, CPG_ISLAND_HG38 n where m.id = n.id and m.id in (select a.id from CPG_ISLAND_HG38 a, REFSEQ_UNIQUE_HG38 b where b.gene_name = '", gene_name, "' ", "and a.chrom = b.chrom and case when b.strand = '+' then ( (a.chromstart >= (b.txstart-3000) and a.chromstart <= b.txend+3000) or (a.chromend >= (b.txstart-3000) and a.chromend <= b.txend+3000) or (a.chromstart <= (b.txstart-3000) and a.chromend >= b.txend+3000) ) else ( (a.chromstart >= b.txstart-3000 and a.chromstart <= b.txend+3000) or (a.chromend >= b.txstart-3000 and a.chromend <= b.txend+3000) or (a.chromstart-3000 <= b.txstart and a.chromend >= b.txend+3000)) end)")
df2 <- dbGetQuery(DB, statement2)

NEPC <- subset(df2, sampleName%in%NEPC_sample)
CRPC <- subset(df2, sampleName%in%CRPC_sample)


statement3 <- paste0("select a.chromstart, a.name, a.chromend, a.chrom, b.strand, b.txstart, b.txend from CPG_ISLAND_HG38 a, REFSEQ_UNIQUE_HG38 b where b.gene_name = '", gene_name, "' ", "and a.chrom = b.chrom and case when b.strand = '+' then ( (a.chromstart >= (b.txstart-3000) and a.chromstart <= b.txend+3000) or (a.chromend >= (b.txstart-3000) and a.chromend <= b.txend+3000) or (a.chromstart <= (b.txstart-3000) and a.chromend >= b.txend+3000) ) else ( (a.chromstart >= b.txstart-3000 and a.chromstart <= b.txend+3000) or (a.chromend >= b.txstart-3000 and a.chromend <= b.txend+3000) or (a.chromstart <= b.txstart-3000 and a.chromend >= b.txend+3000)) end")
df3 <- dbGetQuery(DB, statement3)

dbDisconnect(DB)


# Set the axis labels and ranges
xlim <- range(NEPC$pos, CRPC$pos)
ylim <- c(0, 100)
xlab <- paste0("position on ", df3[1,4])
ylab <- "methylation frac"
tss_pos <- ifelse(df3[1,5]=='+', df3[1,6], df3[1,7])
tts_pos <- ifelse(df3[1,5]=='+', df3[1,7], df3[1,6])

filename = paste0(gene_name, "_beltran_NEPC", "_CRPC", "_cpg_islands")
pdf(file = paste0(filename, "_level.pdf"), onefile=FALSE)

##find the overlapped points
overlapped <- intersect(paste(NEPC$pos, NEPC$frac), paste(CRPC$pos, CRPC$frac))
overlapped_points <- data.frame(matrix(overlapped, ncol = 1, byrow = TRUE), stringsAsFactors = FALSE)
colnames(overlapped_points) = "overlapped"
overlapped_points <- separate(overlapped_points, overlapped, into = c("pos", "frac"), sep = " ")

# Create an empty plot with the appropriate axis labels and ranges
par(cex.axis=0.8)
plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

points(NEPC$pos, NEPC$frac, type="p", col = "skyblue", pch = 16, lwd = 3, cex = 0.5)
points(CRPC$pos, CRPC$frac, type="p", col = "red", pch = 16, lwd = 3, cex = 0.5)
points(overlapped_points$pos, overlapped_points$frac, type="p", col = "darkviolet", pch = 16, lwd = 3, cex = 0.5)


rect(df3[,1], -1, df3[,3], 1.5, col = "lightgreen")
text((df3[,1] + df3[,3]) / 2, 0, labels = c(df3[,2]), cex = 0.3, pos=1)
abline(v = tss_pos, lty = 2, lwd = 2, col="blue")
text(x = tss_pos+500, y = 60, labels = "TSS", col = "blue")
abline(v = tts_pos, lty = 2, lwd = 2, col="black")
text(x = tts_pos+500, y = 60 , labels = "TTS", col = "black")

if(df3[1,5]=='+')
{
	arrows(x0 = tss_pos+200, y0 = 0, x1 = tss_pos+300, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos+400, y0 = 0, x1 = tss_pos+500, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos+600, y0 = 0, x1 = tss_pos+700, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
}else
{
	arrows(x0 = tss_pos-200, y0 = 0, x1 = tss_pos-300, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos-400, y0 = 0, x1 = tss_pos-500, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos-600, y0 = 0, x1 = tss_pos-700, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
}

title(gene_name)
legend("left", legend = c("CRPC", "NEPC", "overlapped"), col = c("red", "skyblue", "darkviolet"), pch = 16, bty = "n")
     
dev.off()



###########################################boxplot cpg_island only

library(patchwork)
plots <- list()
i = 1
for (cpgid in unique(df2[,1])) {
	NEPC_group_all <- NEPC %>% filter(ID==cpgid & disToIsland==0)
	CRPC_group_all <- CRPC %>% filter(ID==cpgid & disToIsland==0)

	NEPC_group <- data.frame(group = "NEPC", value = NEPC_group_all$frac)
	CRPC_group <- data.frame(group = "CRPC", value = CRPC_group_all$frac)
	combined_group <- bind_rows(NEPC_group, CRPC_group)

	my_comparisons <- list(c("NEPC", "CRPC"))

#	ggb = ggviolin(combined_group, x = "group", y = "value", order = c("NEPC", "CRPC"), shape="group", 
#		       color = "black", fill = "group", palette = c("skyblue", "red"), title =  unique(NEPC_group_all$NAME), xlab = "", ylab = "methylation frac") +
#	      geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA)
		       
	ggb = ggboxplot(combined_group, x = "group", y = "value", order = c("NEPC", "CRPC"), add = "jitter", shape="group", 
	       color = "group", palette = c("skyblue", "red"), title =  unique(NEPC_group_all$NAME), xlab = "", ylab = "methylation frac", add.params = list(size = 0.5)) 

	ggb = ggb + stat_compare_means(comparisons = my_comparisons, method = "t.test", label.y = max(combined_group$value)/2) + theme(legend.position = "none")
	plots[[i]] <- ggb
	i = i + 1
}

pdf(file = paste0(filename, "_level_cpgisland_boxplot", ".pdf"), onefile=FALSE)
combined_plot <- wrap_plots(plots, ncol = 3)
combined_plot <- plot_layout(combined_plot + plot_annotation(title = gene_name))
print(combined_plot)
dev.off()














##########################################################################library(RSQLite)
library(dplyr)
library("ggpubr")
library("tidyr")

setwd("D:\\bigdata\\dimelo\\new_code")    

gene_name = 'NKX2-1'

DB <- dbConnect(SQLite(), dbname = "D:\\pdb\\database\\METHYLATION.db")
NEPC_sample <- c('SRR10251320', 'SRR10251334', 'SRR10251336', 'SRR10251374', 'SRR10251401')
CRPC_sample <- c('SRR10251321','SRR10251322','SRR10251323','SRR10251324','SRR10251325','SRR10251326','SRR10251327','SRR10251328','SRR10251329','SRR10251330','SRR10251331','SRR10251332','SRR10251333','SRR10251335','SRR10251337','SRR10251338','SRR10251339','SRR10251340','SRR10251341','SRR10251342','SRR10251344','SRR10251345','SRR10251346','SRR10251347','SRR10251348','SRR10251349','SRR10251350','SRR10251351','SRR10251352','SRR10251353','SRR10251354','SRR10251355','SRR10251356','SRR10251357','SRR10251358','SRR10251359','SRR10251360','SRR10251361','SRR10251362','SRR10251363','SRR10251364','SRR10251365','SRR10251366','SRR10251367','SRR10251368','SRR10251369','SRR10251370','SRR10251371','SRR10251372','SRR10251373','SRR10251375','SRR10251376','SRR10251378','SRR10251379','SRR10251380','SRR10251381','SRR10251382','SRR10251383','SRR10251384','SRR10251385','SRR10251386','SRR10251387','SRR10251388','SRR10251389','SRR10251390','SRR10251391','SRR10251392','SRR10251393','SRR10251394','SRR10251395','SRR10251396','SRR10251397','SRR10251398','SRR10251399','SRR10251400','SRR10251402','SRR10251403','SRR10251404','SRR10251405','SRR10251406','SRR10251407','SRR10251408','SRR10251409','SRR10251410','SRR10251411','SRR10251412','SRR10251413','SRR10251414','SRR10251415','SRR10251416','SRR10251417','SRR10251418','SRR10251419')

statement2 <- paste0("SELECT n.id, n.name, disToIsland, readStart pos, percent frac, sampleName from CPG_ISLAND_FELIX_WGBS m, CPG_ISLAND_HG38 n where m.id = n.id and m.id in (select a.id from CPG_ISLAND_HG38 a, REFSEQ_UNIQUE_HG38 b where b.gene_name = '", gene_name, "' ", "and a.chrom = b.chrom and case when b.strand = '+' then ( (a.chromstart >= (b.txstart-3000) and a.chromstart <= b.txend+3000) or (a.chromend >= (b.txstart-3000) and a.chromend <= b.txend+3000) or (a.chromstart <= (b.txstart-3000) and a.chromend >= b.txend+3000) ) else ( (a.chromstart >= b.txstart-3000 and a.chromstart <= b.txend+3000) or (a.chromend >= b.txstart-3000 and a.chromend <= b.txend+3000) or (a.chromstart-3000 <= b.txstart and a.chromend >= b.txend+3000)) end)")
df2 <- dbGetQuery(DB, statement2)

NEPC <- subset(df2, sampleName%in%NEPC_sample)
CRPC <- subset(df2, sampleName%in%CRPC_sample)


statement3 <- paste0("select a.chromstart, a.name, a.chromend, a.chrom, b.strand, b.txstart, b.txend from CPG_ISLAND_HG38 a, REFSEQ_UNIQUE_HG38 b where b.gene_name = '", gene_name, "' ", "and a.chrom = b.chrom and case when b.strand = '+' then ( (a.chromstart >= (b.txstart-3000) and a.chromstart <= b.txend+3000) or (a.chromend >= (b.txstart-3000) and a.chromend <= b.txend+3000) or (a.chromstart <= (b.txstart-3000) and a.chromend >= b.txend+3000) ) else ( (a.chromstart >= b.txstart-3000 and a.chromstart <= b.txend+3000) or (a.chromend >= b.txstart-3000 and a.chromend <= b.txend+3000) or (a.chromstart <= b.txstart-3000 and a.chromend >= b.txend+3000)) end")
df3 <- dbGetQuery(DB, statement3)

dbDisconnect(DB)


# Set the axis labels and ranges
xlim <- range(NEPC$pos, CRPC$pos)
ylim <- c(0, 100)
xlab <- paste0("position on ", df3[1,4])
ylab <- "methylation frac"
tss_pos <- ifelse(df3[1,5]=='+', df3[1,6], df3[1,7])
tts_pos <- ifelse(df3[1,5]=='+', df3[1,7], df3[1,6])

filename = paste0(gene_name, "_felix_NEPC", "_CRPC", "_cpg_islands")
pdf(file = paste0(filename, "_level.pdf"), onefile=FALSE)

##find the overlapped points
overlapped <- intersect(paste(NEPC$pos, NEPC$frac), paste(CRPC$pos, CRPC$frac))
overlapped_points <- data.frame(matrix(overlapped, ncol = 1, byrow = TRUE), stringsAsFactors = FALSE)
colnames(overlapped_points) = "overlapped"
overlapped_points <- separate(overlapped_points, overlapped, into = c("pos", "frac"), sep = " ")

# Create an empty plot with the appropriate axis labels and ranges
par(cex.axis=0.8)
plot(0, 0, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)

points(NEPC$pos, NEPC$frac, type="p", col = "skyblue", pch = 16, lwd = 3, cex = 0.5)
points(CRPC$pos, CRPC$frac, type="p", col = "red", pch = 16, lwd = 3, cex = 0.5)
points(overlapped_points$pos, overlapped_points$frac, type="p", col = "darkviolet", pch = 16, lwd = 3, cex = 0.5)


rect(df3[,1], -1, df3[,3], 1.5, col = "lightgreen")
text((df3[,1] + df3[,3]) / 2, 0, labels = c(df3[,2]), cex = 0.3, pos=1)
abline(v = tss_pos, lty = 2, lwd = 2, col="blue")
text(x = tss_pos+500, y = 60, labels = "TSS", col = "blue")
abline(v = tts_pos, lty = 2, lwd = 2, col="black")
text(x = tts_pos+500, y = 60 , labels = "TTS", col = "black")

if(df3[1,5]=='+')
{
	arrows(x0 = tss_pos+200, y0 = 0, x1 = tss_pos+300, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos+400, y0 = 0, x1 = tss_pos+500, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos+600, y0 = 0, x1 = tss_pos+700, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
}else
{
	arrows(x0 = tss_pos-200, y0 = 0, x1 = tss_pos-300, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos-400, y0 = 0, x1 = tss_pos-500, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
	arrows(x0 = tss_pos-600, y0 = 0, x1 = tss_pos-700, y1 = 0, length = 0.1, angle = 30, lwd = 2, col="blue")
}

title(gene_name)
legend("left", legend = c("CRPC", "NEPC", "overlapped"), col = c("red", "skyblue", "darkviolet"), pch = 16, bty = "n")
     
dev.off()



###########################################boxplot cpg_island only

library(patchwork)
plots <- list()
i = 1
for (cpgid in unique(df2[,1])) {
	NEPC_group_all <- NEPC %>% filter(ID==cpgid & disToIsland==0)
	CRPC_group_all <- CRPC %>% filter(ID==cpgid & disToIsland==0)

	NEPC_group <- data.frame(group = "NEPC", value = NEPC_group_all$frac)
	CRPC_group <- data.frame(group = "CRPC", value = CRPC_group_all$frac)
	combined_group <- bind_rows(NEPC_group, CRPC_group)

	my_comparisons <- list(c("NEPC", "CRPC"))

#	ggb = ggviolin(combined_group, x = "group", y = "value", order = c("NEPC", "CRPC"), shape="group", 
#		       color = "black", fill = "group", palette = c("skyblue", "red"), title =  unique(NEPC_group_all$NAME), xlab = "", ylab = "methylation frac") +
#	      geom_boxplot(width = 0.05, fill = "white", outlier.shape = NA)
		       
	ggb = ggboxplot(combined_group, x = "group", y = "value", order = c("NEPC", "CRPC"), add = "jitter", shape="group", 
	       color = "group", palette = c("skyblue", "red"), title =  unique(NEPC_group_all$NAME), xlab = "", ylab = "methylation frac", add.params = list(size = 0.5)) 

	ggb = ggb + stat_compare_means(comparisons = my_comparisons, method = "t.test", label.y = max(combined_group$value)/2) + theme(legend.position = "none")
	plots[[i]] <- ggb
	i = i + 1
}

pdf(file = paste0(filename, "_level_cpgisland_boxplot", ".pdf"), onefile=FALSE)
combined_plot <- wrap_plots(plots, ncol = 3)
combined_plot <- plot_layout(combined_plot + plot_annotation(title = gene_name))
print(combined_plot)
dev.off()


