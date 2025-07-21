library(tidyverse)
library(ggplot2)
library(gridExtra)
library(data.table)
library(scales)

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

dt_combinden = setDT(rbind(HA_cov,NA_cov,M_cov))

dt_combinden[, c("a", "b","subtype","d") := tstrsplit(segment, " ", fixed=TRUE)]
dt_combinden[,('segment') := factor(fcase(a == 'H1', 'HA',
                                   a == 'H3', 'HA',
                                   a == 'N1', 'NA',
                                   a == 'N2', 'NA',
                                   a == 'M', 'M'),
                                   levels = c('HA','NA','M'))]
dt_combinden

mybreaks <- function(x) {
  x <- breaks_pretty(n=3)(x)
  x[x==0] <- NA
  x
}

plot_alternative =
  dt_combinden %>%
  ggplot(aes(x=Position,y=Coverage, colour = segment))+
  geom_line(stat = 'identity',size = 1)+
  scale_y_log10("Read Depth",
                breaks = trans_breaks("log10", function(x) 10^x[x != 0], n=3),
                labels = trans_format("log10", math_format(10^.x)))+
  scale_x_continuous("Position", 
                     breaks = mybreaks,
                     guide = guide_axis(angle = 30)
                     )+
  facet_grid(rows = vars(subtype), 
             cols = vars(segment),
             scales = "free_x")+
  theme(legend.title = element_blank(),
        legend.text = element_text(size = 14),
        legend.key.size = unit(1, "lines"),
        legend.position = "",
        axis.title.x = element_text(size = 16,margin = margin(t=15)),
        strip.text.x = element_text(size = 14),
        strip.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, margin = margin(r = 15)),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "lightgrey", 
                                    fill=NA, 
                                    linewidth=1),
        plot.margin = margin(l=10,r=20,t=10,b=10))

plot_alternative



ggsave("/Users/anjohn/Desktop/manuscripts/Influenza/talk/clinical_coverage_alternative.png",
       plot = plot_alternative, width = 7.4, height = 4.4, units = "in", dpi = 600)
