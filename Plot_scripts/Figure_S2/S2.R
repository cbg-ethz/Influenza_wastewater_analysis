library('reader')

# setwd("/Users/anjohn/Desktop/euler/wastewater/Influenza/GitHub/IAV_wastewater/Plot_scripts/Figure_S2/")

#generated using 1) ../../readLength.sh & 2) ../../combining_readLengths.sh
read_len_ex1 = fread('../../rawData_summary/experimet_1_readLengths_combined.tsv.gz')
read_len_ex2 = fread('../../rawData_summary/experimet_2_readLengths_combined.tsv.gz')


dt=bind_rows(read_len_ex1,read_len_ex2)


p <- ggplot(dt, aes(x = V1)) + 
  geom_histogram(color = "blue",fill="blue",binwidth = 25,alpha=0.5) +  # Adjust width to control spacing between bars
  labs(x = "Read Length", y = "Read Count") +
  theme(axis.title.x = element_text(size = 18, margin = margin(t=15)), 
        axis.title.y = element_text(size = 18, margin = margin(r=15)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  coord_cartesian(xlim = c(0,1500))

p

# ggsave(paste0("ReadLen_Distribution.png"), 
#        plot = p, width = 8, height = 4, units = "in", dpi = 600)








#####Clinical samples #####

#generated using ../../readLength.sh
dt_H1 = fread('../../rawData_summary/H1N1_extract_readLengths.txt')
dt_H3 = fread('../../rawData_summary/H3N2_extract_readLengths.txt')


p_H1 <- ggplot(dt_H1, aes(x = V1)) + 
  geom_histogram(color = "blue",fill="blue",binwidth = 25,alpha=0.5) +  # Adjust width to control spacing between bars
  labs(x = "Read Length", y = "Read Count") +
  theme(axis.title.x = element_text(size = 18, margin = margin(t=15)), 
        axis.title.y = element_text(size = 18, margin = margin(r=15)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  coord_cartesian(xlim = c(0,1500))

p_H1

p_H3 <- ggplot(dt_H3, aes(x = V1)) + 
  geom_histogram(color = "blue",fill="blue",binwidth = 25,alpha=0.5) +  # Adjust width to control spacing between bars
  labs(x = "Read Length", y = "Read Count") +
  theme(axis.title.x = element_text(size = 18, margin = margin(t=15)), 
        axis.title.y = element_text(size = 18, margin = margin(r=15)),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16))+
  coord_cartesian(xlim = c(0,1500))

p_H3

