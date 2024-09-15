## each time, just change the name of gene which needs to be checked
## note: need to use R 64 bits, not R 32 bits program

geneA = 'gene1'
geneB = 'gene2'

gene_exp_table = ''

library(scales)
library(RODBC)

library(ReporteRs)
library(ggplot2)
#library(formattable)
library(plyr)
library(htmltools)
library(webshot)    
library(rlist)
library(DT)
library(htmlwidgets)
#library(webshot)
##need install wkhtmltox in c: drive, use this command wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe (or wkhtmltopdf in the future) , don't use webshot since summary table gets bigger, screenshot will not work

## download phantomjs and put this    C:\phantomjs\bin  in the system path variable
## for first summary slide, html table in order to use red color to highlight value

##export_formattable <- function(f, file, width = "100%", height = NULL, 
##                               background = "white", delay = 0.2)
##    {
##      w <- as.htmlwidget(f, width = width, height = height)
##      path <- html_print(w, background = background, viewer = NULL)
##      url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
##      webshot(url,
##              file = file,
##              selector = ".formattable_widget",
##              delay = delay)
##    }
    
writeppt <- function(doc, genecomb, filename, type)
{
	patn = glob2rx(paste0(genecomb, filename))

	file_gene=list.files(path = getwd(), pattern =patn, full.names = TRUE, recursive = TRUE)
	for (file in file_gene)
	{
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, file, width = 7, height = 6.8 , offx = 2, offy = 0.2)
		doc = addFooter(doc, basename(file))
	}
	writeDoc(doc, paste0(genecomb, " correlation ", type, " samples in PCa dataset.pptx"))
}



    
con <- odbcConnect("XXXX", uid="XXXX", pwd="XXXXX", believeNRows=FALSE)
table_list <- sqlQuery(con, "SELECT distinct table_name FROM LAB_DATA_TABLE ORDER BY TABLE_NAME")

summary_table <- data.frame(matrix(ncol = 6, nrow = 0), stringsAsFactors = FALSE)
colnames(summary_table) <- c("dataset", "All", "Normal", "Tumor", "Metastasis", "T_and_M")

geneA = toupper(geneA)
geneB = toupper(geneB)
genecomb = paste0(geneA, "_", geneB)

for (j in 1:length(table_list$TABLE_NAME))
{
	gene_exp_table = table_list$TABLE_NAME[j]
	s1 = paste0("SELECT SAMPLE_ID FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "'")
	s2 = paste0("SELECT SAMPLE_ID FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "' AND TYPE = 'N'")
	s3 = paste0("SELECT SAMPLE_ID FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "' AND TYPE = 'T'")
	s4 = paste0("SELECT SAMPLE_ID FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "' AND TYPE = 'M'")
	s5 = paste0("SELECT SAMPLE_ID FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "' AND TYPE IN ('T', 'M')")
	
	## to check sampleid of s5 are from T+M
	statement_type = paste0("SELECT DISTINCT TYPE FROM LAB_DATA_TABLE WHERE TABLE_NAME = '", gene_exp_table, "' AND TYPE IS NOT NULL")
	ctype <- sqlQuery(con, statement_type)
	stype = as.list(levels(ctype$TYPE))

	datasetname = substr(gene_exp_table, 5, nchar(toString(gene_exp_table)))
	summary_record = list(dataset=datasetname)
	
	for (statement in c(s1, s2, s3, s4, s5))
	{
		sampleid <- sqlQuery(con, statement)

		## add comma
		if (nrow(sampleid))
		{
			pid = paste(shQuote(t(sampleid[,1]), type="sh"), collapse=", ")
			##remove single quote
			pid = gsub("'", '', pid)

			queryA = paste0("SELECT ", pid," FROM ", gene_exp_table, " WHERE UPPER(GENE_NAME) = ", "'", geneA, "'")
			queryB = paste0("SELECT ", pid," FROM ", gene_exp_table, " WHERE UPPER(GENE_NAME) = ", "'", geneB, "'")

			sA <- t(sqlQuery(con, queryA))
			sB <- t(sqlQuery(con, queryB))

			if (statement==s1)
			{
				type = "all samples"
			}else if (statement==s2)
			{
				type = "normal samples"
			}else if (statement==s3)
			{
				type = "tumor samples"
			}else if (statement==s4)
			{
				type = "metastasis samples"
			}else if ((statement==s5)&&((paste0(stype[1],stype[2],stype[3])=="MNT")||(paste0(stype[1],stype[2])=="MT")&&(stype[3]=="NULL")))
			{
				type = "tumor and metastasis samples"
			}
			else
			{
				next
			}

			if (length(sA) && length(sB))
			{
				testd=data.frame(cbind(sA[,1], sB[,1]))
				testd=na.omit(testd)
				if(nrow(testd)>2)
				{
					##cor method: c("pearson", "kendall", "spearman"), pearson is default, can change to test for each
					svalues = cor.test(testd[,1], testd[,2], method = "pearson")
					corr = round(svalues$estimate, digits=3)
					if (!is.na(corr))
					{
						if (svalues$p.value<0.001)
						{
							p = format(svalues$p.value, scientific=TRUE, digits=4)
						}
						else
						{
							p = round(svalues$p.value, digits=3)
						}


						if (type == "all samples")
						{
							summary_record = list.append(summary_record, All=corr)
						}
						if (type == "normal samples")
						{
							summary_record = list.append(summary_record, Normal=corr)
						}
						if (type == "tumor samples")
						{
							summary_record = list.append(summary_record, Tumor=corr)
						}
						if (type == "metastasis samples")
						{
							summary_record = list.append(summary_record, Metastasis=corr)
						}
						if (type == "tumor and metastasis samples")
						{
							summary_record = list.append(summary_record, T_and_M=corr)
						}


						rp = paste0("r = ", corr, ",  ", "p = ", p)
						sp <- ggplot(data=testd, aes(X1,X2))
						sp + geom_point() + geom_point(shape=21, size=3, color="gray", fill="red") + stat_smooth(method=lm, color="blue") + 
						xlab(geneA) + ylab(geneB) + annotate("text", x = min(testd$X1), y = Inf, hjust = 0, vjust = 2, label = rp) + ggtitle(paste0(type, " in ", datasetname)) + 
						theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) + theme(axis.title.x = element_text(size = rel(1.8), angle = 00))+ 
						theme(plot.title = element_text(size=16, hjust = 0.5, face="bold"))

						ggsave(paste0(genecomb, " correlation ", type, " in ", datasetname, ".pdf"), useDingbats=FALSE)
						ggsave(paste0(genecomb, " correlation ", type, " in ", datasetname, ".jpg"))
					}
				}
			}
		}
	}
	summary_table <- rbind.fill(summary_table, as.data.frame(summary_record, stringsAsFactors = FALSE))
}
close(con)


##teset = formattable(summary_table, ALL = color_tile("white", "orange"))
##sumslide = formattable(summary_table, list(
##   All = formatter("span", style = x ~ ifelse((round(as.numeric(x), 4) >= 0.3) || (round(as.numeric(x), 4) <= -0.3), style(color = "red", font.weight = "bold"), NA)),
##   Normal = formatter("span", style = x ~ ifelse((round(as.numeric(x), 4) >= 0.3) || (round(as.numeric(x), 4) <= -0.3), style(color = "red", font.weight = "bold"), NA)),
##   Tumor = formatter("span", style = x ~ ifelse((round(as.numeric(x), 4) >= 0.3) || (round(as.numeric(x), 4) <= -0.3), style(color = "red", font.weight = "bold"), NA)),
##   Metastasis = formatter("span", style = x ~ ifelse((round(as.numeric(x), 4) >= 0.3) || (round(as.numeric(x), 4) <= -0.3), style(color = "red", font.weight = "bold"), NA)),
##   T_and_M = formatter("span", style = x ~ ifelse((round(as.numeric(x), 4) >= 0.3) || (round(as.numeric(x), 4) <= -0.3), style(color = "red", font.weight = "bold"), NA))
##))
##export_formattable(f, sumslidefile)


widget_all = datatable(summary_table, list(dom = 't', ordering=F, pageLength = 30), colnames = c('dataset', 'r value in All', 'r value in Normal', 'r value in Tumor', 'r value in Metastasis', 'r value in Tumor and Metastasis')) %>%
formatStyle('All', backgroundColor = 'skyblue') %>% 
formatStyle(names(summary_table), textAlign = 'center', color = styleInterval(c(-0.3, 0.3)-1e-6, c("red", "black", "green")))
saveWidget(widget_all, "temp_all.html", selfcontained = FALSE)
#sumslidefile_all = "correlation summary all.jpg"
#webshot("temp_all.html", file = sumslidefile_all, cliprect = "viewport")
shell("c:/wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe --javascript-delay 1 temp_all.html correlation_summary_all.jpg")

widget_normal = datatable(summary_table, list(dom = 't', ordering=F, pageLength = 30), colnames = c('dataset', 'r value in All', 'r value in Normal', 'r value in Tumor', 'r value in Metastasis', 'r value in Tumor and Metastasis')) %>%
formatStyle('Normal', backgroundColor = 'skyblue') %>% 
formatStyle(names(summary_table), textAlign = 'center', color = styleInterval(c(-0.3, 0.3)-1e-6, c("red", "black", "green")))
saveWidget(widget_normal, "temp_normal.html", selfcontained = FALSE)
#sumslidefile_normal = "correlation summary normal.jpg"
#webshot("temp_normal.html", file = sumslidefile_normal, cliprect = "viewport")
shell("c:/wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe --javascript-delay 1 temp_normal.html correlation_summary_normal.jpg")

widget_tumor = datatable(summary_table, list(dom = 't', ordering=F, pageLength = 30), colnames = c('dataset', 'r value in All', 'r value in Normal', 'r value in Tumor', 'r value in Metastasis', 'r value in Tumor and Metastasis')) %>%
formatStyle('Tumor', backgroundColor = 'skyblue') %>%
formatStyle(names(summary_table), textAlign = 'center', color = styleInterval(c(-0.3, 0.3)-1e-6, c("red", "black", "green")))
saveWidget(widget_tumor, "temp_tumor.html", selfcontained = FALSE)
#sumslidefile_tumor = "correlation summary tumor.jpg"
#webshot("temp_tumor.html", file = sumslidefile_tumor, cliprect = "viewport")
shell("c:/wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe --javascript-delay 1 temp_tumor.html correlation_summary_tumor.jpg")

widget_metastasis = datatable(summary_table, list(dom = 't', ordering=F, pageLength = 30), colnames = c('dataset', 'r value in All', 'r value in Normal', 'r value in Tumor', 'r value in Metastasis', 'r value in Tumor and Metastasis')) %>%
formatStyle('Metastasis', backgroundColor = 'skyblue') %>% 
formatStyle(names(summary_table), textAlign = 'center', color = styleInterval(c(-0.3, 0.3)-1e-6, c("red", "black", "green")))
saveWidget(widget_metastasis, "temp_metastasis.html", selfcontained = FALSE)
#sumslidefile_metastasis = "correlation summary metastasis.jpg"
#webshot("temp_metastasis.html", file = sumslidefile_metastasis, cliprect = "viewport")
shell("c:/wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe --javascript-delay 1 temp_metastasis.html correlation_summary_metastasis.jpg")

widget_t_and_m = datatable(summary_table, list(dom = 't', ordering=F, pageLength = 30), colnames = c('dataset', 'r value in All', 'r value in Normal', 'r value in Tumor', 'r value in Metastasis', 'r value in Tumor and Metastasis')) %>%
formatStyle('T_and_M', backgroundColor = 'skyblue') %>% 
formatStyle(names(summary_table), textAlign = 'center', color = styleInterval(c(-0.3, 0.3)-1e-6, c("red", "black", "green")))
saveWidget(widget_t_and_m, "temp_t_and_m.html", selfcontained = FALSE)
#sumslidefile_t_and_m = "correlation summary t_and_m.jpg"
#webshot("temp_t_and_m.html", file = sumslidefile_t_and_m, cliprect = "viewport")
shell("c:/wkhtmltox/wkhtmltopdf/bin/wkhtmltoimage.exe --javascript-delay 1 temp_t_and_m.html correlation_summary_t_and_m.jpg")


filename = ""
for (filename in c("*correlation all samples in*.jpg", "*correlation normal samples in*.jpg", "*correlation tumor samples in*.jpg", "*correlation metastasis samples in*.jpg", "*correlation tumor and metastasis samples in*.jpg"))
{
	
	if (filename=="*correlation all samples in*.jpg")
	{
		doc =pptx( )
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, "correlation_summary_all.jpg", width = 8.5, height = 6.8, offx = 2, offy = 0.2)
		writeppt(doc, genecomb, filename, "all")
	}
	else if (filename=="*correlation normal samples in*.jpg")
	{
		doc =pptx( )
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, "correlation_summary_normal.jpg", width = 8.5, height = 6.8, offx = 2, offy = 0.2)
		writeppt(doc, genecomb, filename, "normal")
	}
	else if (filename=="*correlation tumor samples in*.jpg")
	{
		doc =pptx( )
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, "correlation_summary_tumor.jpg", width = 8.5, height = 6.8, offx = 2, offy = 0.2)
		writeppt(doc, genecomb, filename, "tumor")
	}
	else if (filename=="*correlation metastasis samples in*.jpg")
	{
		doc =pptx( )
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, "correlation_summary_metastasis.jpg", width = 8.5, height = 6.8, offx = 2, offy = 0.2)
		writeppt(doc, genecomb, filename, "metastasis")
	}
	else if (filename=="*correlation tumor and metastasis samples in*.jpg")
	{
		doc =pptx( )
		doc = addSlide(doc, "Two Content")
		doc = addImage(doc, "correlation_summary_t_and_m.jpg", width = 8.5, height = 6.8, offx = 2, offy = 0.2)
		writeppt(doc, genecomb, filename, "tumor and metastasis")
	}
}


