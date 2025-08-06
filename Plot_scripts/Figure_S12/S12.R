library(data.table)
library(tidyverse)
library(ggplot2)
library(LaplacesDemon)

# setwd("path_to_Figure_S7/")

#### READ IN DATA ####

dt_Fig3B = fread("dt_Fig3B.tsv")

dt_cov = fread("abundance_estimates_vpipe.csv",
               header = TRUE)[,c("sample","week","year","count")]

dt_map = fread("../sample_mapping/samples_Info_timeseries.csv",
               header = TRUE)

##### DATA FORMATTING ####

dt_map$Barcode = paste0("barcode",str_trim(dt_map$Barcode))
dt_map$Date = as.Date(dt_map$Date, format="%Y_%m_%d")
dt_map$week = strftime(as.Date(dt_map$Date, format="%Y-%m-%d"), 
                            format = "%V")
dt_map$week = as.integer(dt_map$week)

dt_cov = merge(dt_cov,
      dt_map,
      by.x=c("sample","week"),
      by.y=c("Barcode","week")
      )


dt = merge(dt_Fig3B,
           dt_cov)


##### PLOTTING ####

my_breaks = c(25,100, 500, 2500)


p <- ggplot(dt %>% filter(Clade == "6B.1A.5a.2a" & count > 5),
       aes(label = week)) +
      geom_point(aes(logit(relative_abundance_gisaid), 
                     logit(relative_abundance_ww_deconv_lin),
                     colour=count),
                 alpha = 0.7,
                 size=3) +
      scale_colour_gradient(name = "read count", 
                            trans = "log",
                            breaks = my_breaks, 
                            labels = my_breaks) +
      geom_smooth(aes(x = logit(relative_abundance_gisaid),
                      y = logit(relative_abundance_ww_deconv_lin)),
                  method = "lm",na.rm = TRUE,
                  color = "black", alpha = 0.3, size = 0.3)+
      facet_wrap(vars(Clade),scales = "free")+
      geom_text(aes(x = logit(relative_abundance_gisaid), 
                    y = logit(relative_abundance_ww_deconv_lin),
                    label = week,
                    size = 14),
                size = 3, vjust = -1, hjust = -0.1, color = "black")+
      xlab("logit(Abundance estimate clinical data)")+
      ylab("logit(Abundance estimate\n wastewater data)")+
      theme(text = element_text(size = 18),
            # legend.position="none",
            # strip.text = element_text(size = 22),
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
            plot.margin = margin(t = 20, r = 20, b = 15, l = 10))

p

ggsave("AbundanceEstimatesH1_Correlation_5a2a_readcount.png", 
       plot = p, width = 8, height = 6., units = "in", dpi = 600)


##### RESIDUAL POLT #####


lm_models <- dt %>%
  distinct() %>% 
  filter(week != 52) %>%
  group_by(Clade) %>%
  nest() %>%
  mutate(model = map(data, ~ lm(logit(relative_abundance_ww_deconv_lin) ~ logit(relative_abundance_gisaid), data = .)),
         n_resid = map_int(model, ~ length(resid(.))),
         residuals = map(model, resid))

lm_models_augmented <- lm_models %>%
  mutate(data = map2(data, residuals, ~ mutate(.x, residual = .y))) %>%
  select(Clade, data) %>%
  unnest(data)




q <- ggplot(lm_models_augmented %>% filter(Clade == "6B.1A.5a.2a" & count > 5)) +
  geom_point(aes(count, residual),
             size=3) +
  facet_wrap(vars(Clade),scales = "free")+
  xlab("Read count")+
  ylab("Residuals")+
  scale_x_log10() +
  theme(text = element_text(size = 18),
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
        plot.margin = margin(t = 20, r = 20, b = 15, l = 10))


# ggsave("AbundanceEstimatesH1_Residual_readcount.png", 
#        plot = q, width = 8, height = 5., units = "in", dpi = 600)
# 
