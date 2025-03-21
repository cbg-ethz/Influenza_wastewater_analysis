library(tidyverse)
library(ggplot2)
library(data.table)
library(DHARMa)
library(glmmTMB)

#TODO: remove before making public
setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_2")
source("helper_functions.R")

#####LOADING DATA#####

H1_cov_ex1 = "../../work-IA_H1/v-pipe/results_experiment_1/H1_coverages_combined.tsv.gz"
H1_cov_ex2 = "../../work-IA_H1/v-pipe/results_experiment_2/H1_coverages_combined.tsv.gz"

N1_cov_ex1 = "../../work-IA_N1/v-pipe/results_experiment_1/N1_coverages_combined.tsv.gz"
N1_cov_ex2 = "../../work-IA_N1/v-pipe/results_experiment_2/N1_coverages_combined.tsv.gz"


H3_cov_ex1 = "../../work-IA_H3/v-pipe/results_experiment_1/H3_coverages_combined.tsv.gz"
H3_cov_ex2 = "../../work-IA_H3/v-pipe/results_experiment_2/H3_coverages_combined.tsv.gz"

N2_cov_ex1 = "../../work-IA_N2/v-pipe/results_experiment_1/N2_coverages_combined.tsv.gz"
N2_cov_ex2 = "../../work-IA_N2/v-pipe/results_experiment_2/N2_coverages_combined.tsv.gz"


M_cov_ex1 = "../../work-IA_M/v-pipe/results_experiment_1/M_coverages_combined.tsv.gz"
M_cov_ex2 = "../../work-IA_M/v-pipe/results_experiment_2/M_coverages_combined.tsv.gz"

ex1_trans_file = "../sample_mapping/samples_Info_experiment_1.csv"
ex2_trans_file = "../sample_mapping/samples_Info_experiment_2.csv"


dt_H1 = bind_rows(formatting_dt(H1_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(H1_cov_ex2,ex2_trans_file,"_ex2"))
setDT(dt_H1)
dt_H1[, ':='(Segment="HA",
             Subtype="H1N1",
             Segment_Subtype="H1")]


dt_N1 = bind_rows(formatting_dt(N1_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(N1_cov_ex2,ex2_trans_file,"_ex2"))
setDT(dt_N1)
dt_N1[, ':='(Segment="NA",
             Subtype="H1N1",
             Segment_Subtype="N1")]


dt_H3 = bind_rows(formatting_dt(H3_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(H3_cov_ex2,ex2_trans_file,"_ex2"))
setDT(dt_H3)
dt_H3[, ':='(Segment="HA",
             Subtype="H3N2",
             Segment_Subtype="H3")]

dt_N2 = bind_rows(formatting_dt(N2_cov_ex1,ex1_trans_file,"_ex1"),
                  formatting_dt(N2_cov_ex2,ex2_trans_file,"_ex2"))
setDT(dt_N2)
dt_N2[, ':='(Segment="NA",
             Subtype="H3N2",
             Segment_Subtype="N2")]


dt_M = bind_rows(formatting_dt(M_cov_ex1,ex1_trans_file,"_ex1"),
                 formatting_dt(M_cov_ex2,ex2_trans_file,"_ex2"))
setDT(dt_M)
dt_M[, ':='(Segment="M",
            Subtype="M_both",
            Segment_Subtype="none")]



dt_cov = bind_rows(dt_H1,dt_N1,dt_H3,dt_N2,dt_M)

dt = dt_cov %>% 
  group_by(name,Season,Date,Place,concentration,Segment,Subtype,Segment_Subtype,inhibition) %>% 
  summarise(mean_pos_coverag = mean(Coverage)) %>% 
  ungroup() %>% 
  setDT()


##### CLEANING DT ####

dt = dt[, c("Barcode", "Experiment") := tstrsplit(name, "_", fixed = TRUE)]

dt_model = dt
dt_model$name = as.factor(dt_model$name)
dt_model$Subtype = factor(dt_model$Subtype, 
                             levels = c("H1N1",'H3N2', 'M_both'))
dt_model$Segment = factor(dt_model$Segment, 
                             levels = c("HA",'NA', 'M'))
dt_model$Season = factor(dt_model$Season, 
                          levels = c("2022/2023", "2023/2024"))
dt_model$Place = as.factor(dt_model$Place)
dt_model$Experiment = as.factor(dt_model$Experiment)
dt_model$mean_pos_coverag_log = log(dt_model$mean_pos_coverag)
dt_model$mean_pos_coverag_round = round(dt_model$mean_pos_coverag)

dt_model$concentration_log = log(dt_model$concentration)
dt_model[["week"]] = strftime(as.Date(dt_model$Date, format="%Y_%m_%d"), 
                              format = "%V")



dt_model$conc_transformed = dt_model$concentration * 2


##### MODELLING ####

#Zero-inflation visualized

ggplot(dt_model, aes(mean_pos_coverag_round)) +
  geom_histogram() +
  facet_grid(Subtype ~ Segment, margins=TRUE, scales="free")


fit.nb <- glmmTMB(mean_pos_coverag_round ~ concentration + 
                  Subtype +
                  Segment:Subtype+
                  Experiment+
                  Place+
                  (1 | name) +
                  (1 | Season),
                data = dt_model, 
                ziformula = ~ 1,
                family = nbinom2)

summary(fit.nb)


# fit.nb2 <- glmmTMB(mean_pos_coverag_round ~ concentration2 + 
#                     Subtype +
#                     Segment:Subtype+
#                     Experiment+
#                     Place+
#                     (1 | name) +
#                     (1 | Season),
#                   data = dt_model, 
#                   ziformula = ~ 1,
#                   family = nbinom2)
# 
# summary(fit.nb2)






##### TESTING MODEL ASSUMPTIONS #####


simulationOutput <- simulateResiduals(fittedModel = fit.nb, plot = T)

#Dispersion of the fitted model is not significant
testDispersion(fit.nb)
testZeroInflation(simulationOutput)
testResiduals(fit.nb)

#HOMOGENEITY OF VARIANCE IS MET FOR RANDOM EFFECTS
testCategorical(simulationOutput,dt_model$name)
testCategorical(simulationOutput,dt_model$Season)


# simulationOutput <- simulateResiduals(fittedModel = fit.nb, plot = T)




###FIXED EFFECTS

#ex1 = 120, ex2 = 90
# plotResiduals(simulationOutput, dt_model$Experiment)
# 5 observations per sample --> not much, variation expected
plotResiduals(simulationOutput, dt_model$Place)
plotResiduals(simulationOutput, dt_model$Subtype)
plotResiduals(simulationOutput, dt_model$Segment)
plotResiduals(simulationOutput, dt_model$Experiment)

###RANDOM EFFECTS
#2022/2023 = 120, 2023/2024 = 90
plotResiduals(simulationOutput, dt_model$Season)
# plotResiduals(simulationOutput, dt_model$Experiment)
plotResiduals(simulationOutput, dt_model$name)





