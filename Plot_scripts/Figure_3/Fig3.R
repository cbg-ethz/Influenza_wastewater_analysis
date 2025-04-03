library(data.table)
library(tidyverse)
library(boot)
library(ggplot2)
library(RColorBrewer)



setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_3")

#output of DataCleaning.R
dt = fread("dt_test.tsv")

out_dir = "/Users/anjohn/Desktop/manuscripts/Influenza/Figures/"


##### Small amount of cleaning ####

dt[["year"]] = factor(dt[["year"]],
                          levels=c("2023","2024"))

dt[["week"]] = factor(as.character(dt[["week"]]),
                       levels= c("49","50","51","52",
                                 "1","2","3","4",
                                 "5","6","7","8",
                                 "9","10","11","12"))

#Format for proportional area chart

dt_plot_long <- melt(dt, 
                    id.vars = c("week", "Clade", "year"), 
                    measure.vars = c("relative_abundance_ww_deconv_lin", "relative_abundance_gisaid"), 
                    variable.name = "data_source", 
                    value.name = "relative_abundance")

dt_plot_long$data_source <- factor(dt_plot_long$data_source, 
                                  levels = c("relative_abundance_ww_deconv_lin", 
                                             "relative_abundance_gisaid"), 
                                  labels = c("Wastewater Data Geneva", "Clinical Data Europe"))


#hack for coloring the unassigned fraction

dt_unassaigned = dt_plot_long %>% 
  group_by(week,year,data_source) %>% 
  summarise(relative_abundance = round(1-sum(relative_abundance),3)) %>% 
  setDT()

dt_unassaigned[["Clade"]] = "unassigned"

#Wastewater
median_unassigned = median(dt_unassaigned[data_source == 'Wastewater Data Geneva'][['relative_abundance']]) 
mean_unassigned = mean(dt_unassaigned[data_source == 'Wastewater Data Geneva'][['relative_abundance']])



dt_plot = bind_rows(dt_plot_long,dt_unassaigned) %>% 
  filter(week != 52) #drop-out
dt_plot[["Clade"]] = factor(dt_plot[["Clade"]],
                      levels= c("unassigned", "6B.1A.5a.2a","6B.1A.5a.2a.1"))


colour_palette<-c(
  "unassigned"="grey",
  "6B.1A.5a.2a.1"="#1F78B4",
  "6B.1A.5a.2a"="#A6CEE3")





Fig3A <- ggplot(dt_plot, aes(group=Clade))+
  geom_area(aes(x = week, 
                y = relative_abundance,
                fill = Clade),
            colour="white")+
  ylab("Abundance estimates")+
  xlab("Week (2023/2024)")+
  facet_grid(cols=vars(data_source))+
  scale_fill_manual(values = colour_palette,
                    breaks=c("6B.1A.5a.2a",
                             "6B.1A.5a.2a.1",
                             "unassigned"))+
  theme(text = element_text(size = 20),
        legend.title =element_blank(),
        legend.position="bottom",
        # strip.text = element_text(size = 22),
        legend.text =  element_text(margin = margin(r = 15, l=5)),
        axis.title.x = element_text(margin = margin(t = 20)),  
        axis.title.y = element_text(margin = margin(r = 20)),
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        size = 0.5),
        panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "lightgrey", 
                                    fill=NA, 
                                    linewidth=1))
Fig3A


lm_models <- dt %>%
  filter(week != 52) %>% 
  group_by(Clade) %>%
  do(model = lm(logit(relative_abundance_ww_deconv_lin) ~ logit(relative_abundance_gisaid), data = .))

lm_models$adj_r2 <- sapply(lm_models$model, function(m) summary(m)$adj.r.squared)
lm_models$r2 <- sapply(lm_models$model, function(m) summary(m)$r.squared)


Fig3B <- ggplot(dt %>% filter(Clade == "6B.1A.5a.2a" & week != 52)) +
  geom_point(aes(logit(relative_abundance_gisaid), 
                 logit(relative_abundance_ww_deconv_lin),
                 colour = Clade),
             alpha = 0.7,
             size=3) +
  scale_color_brewer(palette = "Dark2")+
  geom_smooth(aes(x = logit(relative_abundance_gisaid),
                  y = logit(relative_abundance_ww_deconv_lin)),
              method = "lm",na.rm = TRUE, 
              color = "black", alpha = 0.3, linewidth = 1)+
  facet_wrap(vars(Clade),scales = "free")+
  geom_text(data = lm_models %>% filter(Clade == "6B.1A.5a.2a"),
            aes(x = Inf, y = Inf,
                label = paste("adj. R² = ", round(adj_r2, 2))),
            inherit.aes = FALSE,
            hjust = 1.1, vjust = 1.5, size = 6)+
  xlab("logit(Abundance estimate clinical data)")+
  ylab("logit(Abundance estimate\n wastewater data)")+
  theme(text = element_text(size = 18),
        legend.position="none",
        panel.background = element_rect(fill = "white",
                                        colour = "white",
                                        linewidth = 0.5),
        panel.grid.major = element_line(linewidth = 0.25, linetype = 'solid',
                                        colour = "lightgrey"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "lightgrey", 
                                    fill=NA, 
                                    linewidth=1),
        axis.title.x = element_text(margin = margin(t = 20)),  
        axis.title.y = element_text(margin = margin(r = 20)),
        plot.margin = margin(t = 45, r = 20, b = 15, l = 10))

Fig3B







# ggsave(paste0(out_dir,"FIG3A.png"), 
#        plot = Fig3A, width = 10, height = 5, units = "in", dpi = 600)
# 
# 
# ggsave(paste0(out_dir,"FIG3B.png"), 
#        plot = Fig3B, width = 8, height = 5, units = "in", dpi = 600)





