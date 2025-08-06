library(data.table)
library(tidyverse)

H1_ampCov <- fread("H1_comb_ampCoverage.tsv")
H1_ampCov[, ':='(Segment="HA",
                 Subtype="H1N1",
                 Segment_Subtype="H1")]



H3_ampCov <- fread("H3_comb_ampCoverage.tsv")
H3_ampCov[, ':='(Segment="HA",
                 Subtype="H3N2",
                 Segment_Subtype="H3")]


N1_ampCov <- fread("N1_comb_ampCoverage.tsv")
N1_ampCov[, ':='(Segment="NA",
                 Subtype="H1N1",
                 Segment_Subtype="N1")]


N2_ampCov <- fread("N2_comb_ampCoverage.tsv")
N2_ampCov[, ':='(Segment="NA",
                 Subtype="H3N2",
                 Segment_Subtype="N2")]

HA_ampCov <- rbind(H1_ampCov,H3_ampCov)
NA_ampCov <- rbind(N1_ampCov,N2_ampCov)


M_ampCov <- fread("M_comb_ampCoverage.tsv")
M_ampCov[, ':='(Segment="M",
                Subtype="M_both",
                Segment_Subtype="none")]


dt_amp_model = rbind(HA_ampCov,NA_ampCov,M_ampCov)




wilcoxPerAmp <- function(dt_amp){
  
  dt_amp = as.data.table(spread(dt_amp[,c("amplicon",
                                             "read_count",
                                             "segment",
                                             "Date_Location",
                                             "Season")], 
                                 key = amplicon, 
                                 value = read_count))
  
  
  amp = colnames(dt_amp)[str_starts(colnames(dt_amp), 'amp')]
  
  
  results = data.frame(Baselime_amp = character(),
                       amp = character(),
                       V = numeric(),
                       p_value = numeric())
  
  for (i in amp) {
    
    for (j in amp){
      
      if (i == j){
        
        results <- rbind(results, data.frame(
          Baselime_amp = i,
          amp = j,
          V = NA,
          p_value = NA))
        
      }else{
        
        print(i)
        test = wilcox.test(dt_amp[[i]], dt_amp[[j]],
                           paired = TRUE)
        
        results <- rbind(results, data.frame(
          Baselime_amp = i,
          amp = j,
          V = as.numeric(test$statistic),
          p_value = test$p.value
        ))
        
        
      }
    }
  }
  
  results$p_bonferroni <- p.adjust(results$p_value, method = "bonferroni")
  
  results_clean <- results %>%
    na.omit() %>%
    mutate(significance = case_when(
      p_bonferroni < 0.001 ~ "***",
      p_bonferroni < 0.01  ~ "**",
      p_bonferroni < 0.05  ~ "*",
      TRUE                 ~ ""
    ))
  
  
  return(results_clean)
  
  
}



out_table = dt_amp_model %>% 
              mutate(Segment_Subtype = case_when(
                Segment_Subtype == "none" ~ "M",
                .default = Segment_Subtype
              )) %>% 
              group_by(Segment_Subtype,Season)%>%
              do(data.frame(val=wilcoxPerAmp(.))) 


write_csv(out_table,file = "amp_readCount_testing.csv" )
