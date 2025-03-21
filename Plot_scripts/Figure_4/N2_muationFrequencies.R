library(tidyverse)
library(data.table)
library(ggplot2)
library(gridExtra)
library(Biostrings)
library(ggforce)

#CODE DOES NOT WORK AS INTENDED


source("/Users/anjohn/projects/wastewater/influenza/mutationTranslation.R")
source("/Users/anjohn/projects/wastewater/influenza/readingFunctions.R")

N1_resistance_mut = fread("N1_resistanceMuts.tsv")$AA_mutations


dir_euler <- '/Users/anjohn/Desktop/euler/wastewater/Influenza/work-IA_N1/v-pipe/'
segment <- "N1"
fasta_file <- paste0(dir_euler,"references/nexstrain/N1.fasta")
results_aire <-  "results_Aire/"

  
start = 9	
end = 1418
  


samples_aire <- read.table('../sample_mapping/samples_Info_timeseries.csv',
                           sep = ",", header = TRUE)
samples_dir <- paste0("barcode",str_trim(samples_aire$Barcode),"/","batch1")
var_files_aire <- paste0(dir_euler,results_aire,samples_dir,"/variants/SNVs/snvs.vcf")

var_list_aire <- var_files_aire %>% 
  map(function(x) read_and_mark_vcf(x))

var_list_aire <- var_list_aire[!sapply(var_list_aire, is.null)]





NA_mutAA_aire <- var_list_aire %>% 
  map(function(x) identify_AA(fasta_file,start,end,x)) %>% 
  bind_rows() %>%
  merge(.,samples_aire[c('Barcode','Season','Place','Date')],
        by.x=c('barcode'),
        by.y=c('Barcode'))%>% 
  mutate(run = 'rerun')


setDT(NA_mutAA_aire)



# dt <- NA_mutAA_comb[new_AA != org_AA]
# dt[, mut_lable := paste0(org_AA,gene_codonIDX,new_AA)]
# dt[, Date := as.Date(Date, format="%Y_%m_%d")]
# dt[, loc_lable := sapply(str_split(Place, pattern = " "), `[`, 2)]
# 
# 
# p <- ggplot(dt, aes(x = as.factor(mut_lable), y = as.factor(Date), fill = AF)) +
#   geom_tile()+
#   facet_grid(vars(loc_lable), switch = "y", scales = "free") +
#   scale_fill_gradient(low = "lightblue", high = "darkblue")+
#   force_panelsizes(rows = c(0.75,0.22,0.22,0.75))+
#   # scale_y_discrete(breaks=dt_plot$Date, labels=dt_plot$week)+
#   # geom_text(aes(label = round(frac, 2)), color = "white", size = 4) +
#   theme_light()+
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_text(size = 12),
#         axis.text.x = element_text(size = 12,angle = 90, hjust = 1),
#         legend.title = element_text(vjust = .8, size = 12),
#         legend.text = element_text(size = 10),
#         legend.position = "bottom",
#         legend.key.width = unit(1, "cm"),
#         strip.text = element_text(size = 12))+
#   labs(fill = "Mutation Frequency")
# 
# 
# p
# 
# ggsave(paste0("Figures/",segment,"_muatFrequencies.png"), 
#        plot = p, width = 12, height = 8, units = "in", dpi = 300)
# 






#don't care about synonymous mutations
dt_aire <- NA_mutAA_aire[new_AA != org_AA]
dt_aire[, mut_lable := paste0(org_AA,gene_codonIDX,new_AA)]
dt_aire[, Date := as.Date(Date, format="%Y_%m_%d")]
dt_aire[, loc_lable := sapply(str_split(Place, pattern = " "), `[`, 2)]
dt_aire[["week"]] = strftime(as.Date(dt_aire$Date, format="%Y_%m_%d"), 
                             format = "%V")
dt_aire[["year"]] = strftime(as.Date(dt_aire$Date, format="%Y_%m_%d"), 
                             format = "%Y")

dt_aire[["year"]] = factor(dt_aire[["year"]],levels=c("2024","2023"))



q <- ggplot(dt_aire, aes(x = as.factor(gene_codonIDX), y = as.factor(Date), fill = AF)) +
  geom_tile()+
  facet_grid(vars(year), switch = "y", scales = "free") +
  scale_fill_gradient(low = "lightblue", high = "darkblue")+
  force_panelsizes(rows = c(0.8,0.35))+
  scale_x_discrete(breaks=dt_aire$gene_codonIDX, labels=dt_aire$mut_lable)+
  scale_y_discrete(name = "Week",
                   labels = setNames(dt_aire$week, dt_aire$Date))+
  # geom_text(aes(label = round(frac, 2)), color = "white", size = 4) +
  theme_light()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18,margin = margin(r=10)),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16,angle = 90, hjust = 1,vjust = 0.3),
        legend.title = element_text(vjust = .8, size = 16),
        legend.text = element_text(size = 14),
        legend.position = "bottom",
        legend.key.width = unit(1, "cm"),
        strip.text = element_text(size = 16))+
  labs(fill = "Mutation Frequency")

q

ggsave(paste0("Figures/Intermediate_DTs/N2_Aire_muatFrequencies.png"), 
       plot = q, width = 12, height = 8, units = "in", dpi = 300)

