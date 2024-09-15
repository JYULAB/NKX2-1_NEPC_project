library(ggplot2)
library(dplyr)
library(ggpubr)

# Read the data
setwd("D:\\bigdata\\dimelo\\4w_by_guppy")
data<-read.table("FOXA2_4W_PEAKONLY_q4_IN_ALL_75PROB_A_CG_counts.txt", header=T)

# Group by DAY and mA_counts, then summarize to find max mCG_counts
grouped_data <- data %>%
  group_by(DAY, mA_counts) %>%
  summarize(max_mCG_counts = max(mCG_counts), .groups = 'drop')
  
grouped_data$DAY <- as.factor(grouped_data$DAY)

ggscatter(data, x = "mA_counts", y = "mCG_counts", xlim = c(0,80), ylim = c(0,40),
   color = "DAY", palette = "jco",
   ##add = "reg.line", conf.int = TRUE, 
   legend="right", size = 3, alpha = 0.2, position = position_jitter(width = 0.5, height = 0.5)) + 
   guides(color = guide_legend(override.aes = list(size = 6))) + 
   theme(legend.box = "vertical") + 
   geom_smooth(data=grouped_data, aes(x = mA_counts, y = max_mCG_counts, group = DAY, color = DAY),
   method = 'loess', se = FALSE, span= 0.8, linewidth = 1.5) + 
   scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set x-axis to start at zero
   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set y-axis to start at zero
   labs(title = "",
       x = "mA_counts",
       y = "mCG_counts") +
   theme_minimal() +
   theme(legend.title = element_blank()) 
   
   




library(ggplot2)
library(dplyr)
library(ggpubr)

# Read the data
setwd("D:\\bigdata\\dimelo\\nanonome\\LuCaP145_1_DiMeLo")
data<-read.table("LuCaP145_1_peak_q4_q1_75PROB_A_CG_counts_without_outliers.txt", header=T)

# Group by Quantile and mA_counts, then summarize to find max mCG_counts
grouped_data <- data %>%
  group_by(Quantile, mA_counts) %>%
  summarize(max_mCG_counts = max(mCG_counts), .groups = 'drop')
  
grouped_data$Quantile <- as.factor(grouped_data$Quantile)

ggscatter(data, x = "mA_counts", y = "mCG_counts", xlim = c(0,30), ylim = c(0,150),
   color = "Quantile", palette = "jco",
   ##add = "reg.line", conf.int = TRUE, 
   legend="right", size = 3, alpha = 0.2, position = position_jitter(width = 0.5, height = 0.5)) + 
   guides(color = guide_legend(override.aes = list(size = 6))) + 
   theme(legend.box = "vertical") + 
   geom_smooth(data=grouped_data, aes(x = mA_counts, y = max_mCG_counts, group = Quantile, color = Quantile),
   method = 'loess', se = FALSE, span= 0.8, linewidth = 1.5) + 
   scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set x-axis to start at zero
   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set y-axis to start at zero
   labs(title = "",
       x = "mA_counts",
       y = "mCG_counts") +
   theme_minimal() +
   theme(legend.title = element_blank()) 





library(ggplot2)
library(dplyr)
library(ggpubr)

# Read the data
setwd("D:\\bigdata\\dimelo\\nanonome\\LuCaP145_2_DiMeLo")
data<-read.table("LuCaP145_2_peak_q4_q1_75PROB_A_CG_counts_without_outliers.txt", header=T)

# Group by Quantile and mA_counts, then summarize to find max mCG_counts
grouped_data <- data %>%
  group_by(Quantile, mA_counts) %>%
  summarize(max_mCG_counts = max(mCG_counts), .groups = 'drop')
  
grouped_data$Quantile <- as.factor(grouped_data$Quantile)

ggscatter(data, x = "mA_counts", y = "mCG_counts", xlim = c(0,30), ylim = c(0,80),
   color = "Quantile", palette = "jco",
   ##add = "reg.line", conf.int = TRUE, 
   legend="right", size = 3, alpha = 0.2, position = position_jitter(width = 0.5, height = 0.5)) + 
   guides(color = guide_legend(override.aes = list(size = 6))) + 
   theme(legend.box = "vertical") + 
   geom_smooth(data=grouped_data, aes(x = mA_counts, y = max_mCG_counts, group = Quantile, color = Quantile),
   method = 'loess', se = FALSE, span= 0.8, linewidth = 1.5) + 
   scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set x-axis to start at zero
   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set y-axis to start at zero
   labs(title = "",
       x = "mA_counts",
       y = "mCG_counts") +
   theme_minimal() +
   theme(legend.title = element_blank()) 







library(ggplot2)
library(dplyr)
library(ggpubr)

# Read the data
setwd("D:\\bigdata\\dimelo\\nanonome\\NCI_H660_DiMeLo")
data<-read.table("NCI_H660_m1043_peak_q4_q1_75PROB_A_CG_counts_without_outliers.txt", header=T)

# Group by Quantile and mA_counts, then summarize to find max mCG_counts
grouped_data <- data %>%
  group_by(Quantile, mA_counts) %>%
  summarize(max_mCG_counts = max(mCG_counts), .groups = 'drop')
  
grouped_data$Quantile <- as.factor(grouped_data$Quantile)

ggscatter(data, x = "mA_counts", y = "mCG_counts", xlim = c(0,40), ylim = c(0,120),
   color = "Quantile", palette = "jco",
   ##add = "reg.line", conf.int = TRUE, 
   legend="right", size = 3, alpha = 0.2, position = position_jitter(width = 0.5, height = 0.5)) + 
   guides(color = guide_legend(override.aes = list(size = 6))) + 
   theme(legend.box = "vertical") + 
   geom_smooth(data=grouped_data, aes(x = mA_counts, y = max_mCG_counts, group = Quantile, color = Quantile),
   method = 'loess', se = FALSE, span= 0.8, linewidth = 1.5) + 
   scale_x_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set x-axis to start at zero
   scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +  # Set y-axis to start at zero
   labs(title = "",
       x = "mA_counts",
       y = "mCG_counts") +
   theme_minimal() +
   theme(legend.title = element_blank()) 


