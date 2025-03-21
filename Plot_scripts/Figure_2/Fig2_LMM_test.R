library(tidyverse)
library(ggplot2)
library(gridExtra)
library(data.table)
library(lme4)
library(lmerTest)
library(robustlmm)

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
            Subtype="none",
            Segment_Subtype="none")]



dt_cov = bind_rows(dt_H1,dt_N1,dt_H3,dt_N2,dt_M)







dt = dt_cov %>% 
      group_by(name,Season,Date,Place,concentration,Segment,Subtype,Segment_Subtype) %>% 
      summarise(mean_pos_coverag = mean(Coverage)) %>% 
  ungroup() %>% 
  setDT()


dt = dt[, c("Barcode", "Experiment") := tstrsplit(name, "_", fixed = TRUE)]

dt %>% 
  ggplot(aes(x=concentration,y=mean_pos_coverag,color =Subtype))+
  geom_point()+
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)

dt %>% 
  ggplot(aes(x=concentration,y=mean_pos_coverag,color =Place))+
  geom_point()+
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)

dt %>% 
  ggplot(aes(x=concentration,y=mean_pos_coverag,color =Experiment))+
  geom_point()+
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)

dt %>% 
  ggplot(aes(x=concentration,y=mean_pos_coverag,color = Season))+
  geom_point()+
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)




#####TUTORIAL#####
# https://people.math.ethz.ch/~meier/teaching/anova/random-and-mixed-effects-models.html
# https://bookdown.org/steve_midway/DAR/random-effects.html

#name = sample identifier, expect similar behaviour within observations of one 
#sample
# dt[mean_pos_coverag == 0, mean_pos_coverag := 0.001]


stripchart(log(mean_pos_coverag) ~ name, vertical = TRUE, pch = 1, xlab = "sample", 
           data = dt)


dt %>% 
  ggplot(aes(x=name,
             y=mean_pos_coverag,
             col= Subtype,
             group= Subtype))+
  geom_point()+
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  geom_smooth(method = lm,
              se = FALSE)+
  theme(axis.text.x = element_text(angle = 90))


###random intercept model ###
fit.samples <- lmer(mean_pos_coverag ~ (1 | name), data = dt)

summary(fit.samples)

icc = (1212777/(1212777+28950165))*100
#4% of total variance of mean coverage is due to sample difference


confint(fit.samples, oldNames = FALSE)

#checking model assumptions
plot(fit.samples)
par(mfrow = c(1, 2))
qqnorm(ranef(fit.samples)$name[,"(Intercept)"], 
       main = "Random effects")
qqnorm(resid(fit.samples), main = "Residuals")



##2-factor model ####

str(dt)

xtabs(~ name + Subtype, data = dt)
par(mfrow = c(1, 1))
with(dt, interaction.plot(x.factor = name, 
                               trace.factor = Subtype, 
                               response = log(mean_pos_coverag)))

fit.two <- lmer(mean_pos_coverag ~ (1 | Subtype) + (1 | name) + 
                      (1 | Subtype:name), data = dt)

summary(fit.two)

Subtype_var=11420149
name_var=4026379
interaction_var=36466840
error_term=49813


total_var = Subtype_var+name_var+interaction_var+error_term

max(Subtype_var,name_var,interaction_var)
#inconsistency of Subtype coverage across samples (interaction term)
#makes up biggest part of variance (70%)
  
  
confint(fit.two, oldNames = FALSE)
  
  
  
### NESTED MODEL ####
#Place (location of wwtp) could be a batch of some sort

ggplot(dt, aes(y = Subtype , x = mean_pos_coverag)) + 
  geom_point() + 
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)+
  facet_grid(Place ~ .)

fit.nested <- lmer(mean_pos_coverag ~ (1 | Place/Subtype), data = dt)
summary(fit.nested)


### MIXED EFFECT MODEL ####
#for each sample measurements were taken multiple times
#but maybe viral load is more predictive of sample

dt_model = dt
# dt_model$concentration = as.factor(dt_model$concentration)
dt_model$name = as.factor(dt_model$name)
dt_model$Subtype = as.factor(dt_model$Subtype)
dt_model$Place = as.factor(dt_model$Place)
dt_model$mean_pos_coverag_log = log(dt_model$mean_pos_coverag)
dt_model$concentration_log = log(dt_model$concentration)
dt_model[["week"]] = strftime(as.Date(dt_model$Date, format="%Y_%m_%d"), 
                              format = "%V")


ggplot(dt_model, 
       aes(x = concentration, 
               y = mean_pos_coverag, 
               group = Subtype, 
               col = Subtype)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)
  # stat_summary(fun = mean, geom = "line") 

ggplot(dt_model, 
       aes(x = concentration, 
           y = mean_pos_coverag, 
           group = Segment_Subtype, 
           col = Segment_Subtype)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)



ggplot(dt_model, 
       aes(x = concentration, 
           y = mean_pos_coverag, 
           group = Place, 
           col = Place)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)


ggplot(dt_model, 
       aes(x = concentration, 
           y = mean_pos_coverag, 
           group = Season, 
           col = Season)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)

ggplot(dt_model, 
       aes(x = concentration, 
           y = mean_pos_coverag, 
           group = Experiment, 
           col = Experiment)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)

ggplot(dt_model, 
       aes(x = concentration, 
           y = mean_pos_coverag, 
           group = week, 
           col = week)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)



ggplot(dt_model[mean_pos_coverag >= 1, ], 
       aes(x = concentration, 
               y = mean_pos_coverag, 
               group = Segment_Subtype, 
               col = Segment_Subtype)) + 
  geom_point() + 
  scale_y_continuous(trans='log10',
                     labels = scales::scientific)+
  scale_x_continuous(trans='log10',
                     labels = scales::scientific)+
  stat_summary(fun = mean, geom = "line") 

ggplot(dt_model, aes(y = mean_pos_coverag_log, x = concentration_log)) +
  geom_point() +
  xlab("Log(mean coverage)") +
  ylab("Log(concentration)") +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", formula = 'y ~ x', se=F,fullrange = T) +
  facet_wrap(Place~Segment_Subtype)












fit.virLoads <- lmerTest::lmer(log(mean_pos_coverag) ~ concentration + 
                                 (1 | Place) +
                                 (1 | Subtype) +
                       (1 | Subtype:concentration:Place), 
                       data = dt_model)


fit.randSlopes <- lmerTest::lmer(log(mean_pos_coverag) ~ concentration + 
                                 (Subtype | 1) +
                                 (1 | Subtype:concentration:Place), 
                               data = dt_model)


#random intercept
#random effect of Segment_Subtype is nested in Subtype
ggplot(dt_model, aes(y = mean_pos_coverag_log, x = concentration_log)) +
  geom_point() +
  xlab("Log(concentration)") +
  ylab("Log(mean coverage)") +
  theme_classic(base_size = 15) +
  stat_smooth(method = "lm", formula = 'y ~ x', se=F,fullrange = T) +
  facet_wrap(Subtype~Segment)



fit.nested <- lmerTest::lmer(mean_pos_coverag_log ~ concentration_log + 
                            # (concentration_log|Subtype),
                            (1 | Subtype/Segment),
                          data = dt_model)


fit.lmm <-  lmerTest::lmer(mean_pos_coverag_log ~ concentration_log + 
                             Subtype +
                             Subtype:Segment+
                             (1 |Experiment+name+Place+Season),
                           data = dt_model)


rob_fit.lmm <-  rlmer(mean_pos_coverag_log ~ concentration_log + 
                             Subtype +
                             Subtype:Segment+
                             (1 |Experiment+name+Place+Season),
                           data = dt_model)





summary(fit.lmm)
summary(rob_fit.lmm)

lmerTest::anova(rob_fit.lmm)

plot.rlmerMod(rob_fit.lmm)

plot(rob_fit.lmm)


anova(rob_fit.lmm)


par(mfrow = c(2, 3))
plot(fit.lmm)
qqnorm(ranef(fit.lmm)$name[, 1], 
       main = "Random effects of Sample")
qqnorm(ranef(fit.lmm)$Place[, 1], 
       main = "Random effects of Place")
qqnorm(ranef(fit.lmm)$Season[, 1], 
       main = "Random effects of Season")
qqnorm(ranef(fit.lmm)$Experiment[, 1], 
       main = "Random effects of Experiment")
qqnorm(resid(fit.lmm), main = "Residuals")




qqnorm(ranef(fit.lmm)$'Segment:Subtype'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit.lmm), main = "Residuals")





summary(fit.nested)
anova(fit.nested)



## PLOTS fit.nested
plot(fit.nested)
par(mfrow = c(1, 3))
qqnorm(ranef(fit.nested)$Subtype[, 1], 
       main = "Random effects of Subtype")
qqnorm(ranef(fit.nested)$'Segment:Subtype'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit.nested), main = "Residuals")




## PLOTS fit.nested
plot(fit.nested)
par(mfrow = c(1, 3))
qqnorm(ranef(fit.nested)$Subtype[, 1], 
       main = "Random effects of Subtype")
qqnorm(ranef(fit.nested)$'Segment:Subtype'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit.nested), main = "Residuals")



summary(fit.lmm)
anova(fit.lmm)


summary(fit.virLoads)
summary(fit.randSlopes) 
anova(fit.randSlopes)





## PLOTS fit.virLoads
plot(fit.virLoads)
par(mfrow = c(1, 3))
qqnorm(ranef(fit.virLoads)$Subtype[, 1], 
       main = "Random effects of Subtype")
qqnorm(ranef(fit.virLoads)$'Subtype:concentration'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit.virLoads), main = "Residuals")



plot(fit.randSlopes)
par(mfrow = c(1, 3))
qqnorm(ranef(fit.randSlopes)$Subtype[, 1], 
       main = "Random effects of Subtype")
qqnorm(ranef(fit.randSlopes)$'Subtype:concentration'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit.randSlopes), main = "Residuals")










# dt_model_ex2 = dt_model[Experiment=="ex2", ]
# dt_model_ex2$res = resid(fit.virLoads)
# dt_model_ex2$fit = fitted(fit.virLoads)
# dt_model_ex2[res < -1 & fit < -1 , ]




dt_model$res = resid(fit.virLoads)
dt_model$fit = fitted(fit.virLoads)
dt_model[res < -1 & fit < -1 , ]



# identify(fitted_values, residuals, labels = dt_model$name)



# ## PLOTS fit
plot(fit)
par(mfrow = c(1, 3))
qqnorm(ranef(fit)$Subtype[, 1], 
       main = "Random effects of Subtype")
qqnorm(ranef(fit)$'Subtype:concentration'[, 1], 
       main = "Random interaction")
qqnorm(resid(fit), main = "Residuals")











### ZERO-INFLATION MODEL ####
#https://stats.oarc.ucla.edu/r/dae/zinb/
require(pscl)
require(MASS)
require(boot)

dt_model = dt
# dt_model$concentration = as.factor(dt_model$concentration)
dt_model$name = as.factor(dt_model$name)
dt_model$Subtype = as.factor(dt_model$Subtype)
dt_model$Place = as.factor(dt_model$Place)
dt_model$Experiment = as.factor(dt_model$Experiment)
dt_model$Season = as.factor(dt_model$Season)


dt_model$mean_pos_coverag_log = log(dt_model$mean_pos_coverag)
dt_model$concentration_log = log(dt_model$concentration)
dt_model[["week"]] = strftime(as.Date(dt_model$Date, format="%Y_%m_%d"), 
                              format = "%V")

ggplot(dt_model, aes(mean_pos_coverag, fill = Subtype)) +
  geom_histogram() +
  scale_x_log10() +
  facet_grid(Subtype ~ ., margins=TRUE, scales="free_y")



m1 <- zeroinfl(mean_pos_coverag ~ concentration + 
                 Subtype +
                 Subtype:Segment +
                 (1 | Experiment) +
                 (1 | name) +
                 (1 | Place) +
                 (1 | Season),
               data = dt_model, dist = "negbin")
summary(m1)



# install.packages("glmmTMB")
library(glmmTMB)

m1 <- glmmTMB(mean_pos_coverag ~ concentration + 
                Subtype + 
                Subtype:Segment +
                (1 | Experiment) +
                (1 | name) +
                (1 | Place) +
                (1 | Season),
              data = dt_model, 
              family = nbinom2)

summary(m1)
