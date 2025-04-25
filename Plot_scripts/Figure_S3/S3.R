library(tidyverse)
library(ggplot2)
library(gridExtra)

# setwd("path_to_Figure_S3/")
source("helper_functions.R")


##### CLINICAL SAMPLES #####


file_cov_H1 <- '../../work-IA_H1/v-pipe/results_clinical/barcode01/batch1/alignments/coverage.tsv.gz'
file_cov_H3 <- '../../work-IA_H3/v-pipe/results_clinical/barcode01/batch1/alignments/coverage.tsv.gz'

file_cov_N1 <- '../../work-IA_N1/v-pipe/results_clinical/barcode01/batch1/alignments/coverage.tsv.gz'
file_cov_N2 <- '../../work-IA_N2/v-pipe/results_clinical/barcode01/batch1/alignments/coverage.tsv.gz'

file_cov_M_H1N1 <- '../../work-IA_M/v-pipe/results_clinical_H1N1/barcode01/batch1/alignments/coverage.tsv.gz'
file_cov_M_H3N2 <- '../../work-IA_M/v-pipe/results_clinical_H3N2/barcode01/batch1/alignments/coverage.tsv.gz'


cov_H1 <- read_and_mark_cov(file_cov_H1) %>% mutate(segment = "H1 in H1N1 extract")
cov_H3 <- read_and_mark_cov(file_cov_H3) %>% mutate(segment = "H3 in H3N2 extract")

cov_N1 <- read_and_mark_cov(file_cov_N1) %>% mutate(segment = "N1 in H1N1 extract")
cov_N2 <- read_and_mark_cov(file_cov_N2) %>% mutate(segment = "N2 in H3N2 extract")

cov_M_H1 <- read_and_mark_cov(file_cov_M_H1N1) %>% mutate(segment = "M in H1N1 extract")
cov_M_H3 <- read_and_mark_cov(file_cov_M_H3N2) %>% mutate(segment = "M in H3N2 extract")

HA_cov <- rbind(cov_H1,cov_H3)
NA_cov <- rbind(cov_N1,cov_N2)
M_cov <- rbind(cov_M_H1,cov_M_H3)


HA <- HA_cov %>%
  ggplot(aes(x=Position,y=Coverage,colour = segment))+
  geom_line(stat = 'identity',linewidth = 1.5)+
  scale_y_continuous('Read Depth',
                     trans='log10')+
  scale_x_continuous("Position")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(2, "lines"),
        legend.position = "bottom",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.text.x = element_text(size = 18, margin = margin(t = 15)),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(l=10,r=20,t=10,b=10))

# ggsave(paste0("HA_clinical_coverage.png"), 
#        plot = HA, width = 8, height = 4, units = "in", dpi = 600)


Na <- NA_cov %>%
  ggplot(aes(x=Position,y=Coverage,colour = segment))+
  geom_line(stat = 'identity',linewidth = 1.5)+
  scale_y_continuous('Read Depth',
                     trans='log10')+
  scale_x_continuous("Position")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "bottom",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.text.x = element_text(size = 18, margin = margin(t = 15)),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(l=10,r=20,t=10,b=10))

# ggsave(paste0("NA_clinical_coverage.png"), 
#        plot = Na, width = 8, height = 4, units = "in", dpi = 600)



M <- M_cov %>%
  ggplot(aes(x=Position,y=Coverage,colour = segment))+
  geom_line(stat = 'identity',size = 1.5)+
  scale_y_continuous('Read Depth',
                     trans='log10')+
  scale_x_continuous("Position")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 16),
        legend.key.size = unit(1.5, "lines"),
        legend.position = "bottom",
        axis.title.x = element_text(size = 20),
        axis.title.y = element_text(size = 20, margin = margin(r = 15)),
        axis.text.x = element_text(size = 18, margin = margin(t = 15)),
        axis.text.y = element_text(size = 18),
        plot.margin = margin(l=10,r=20,t=10,b=10))

# ggsave(paste0("M_clinical_coverage.png"), 
#        plot = M, width = 8, height = 4, units = "in", dpi = 600)


