library(stringr)
library(data.table)
library(ggplot2)
library(gridExtra)
library(scales)

formatting_dt <- function(file,
                          trans_file,
                          extension=""){
  
  dt = fread(file, col.names = c("Ref","Position","Coverage","sample"))
  
  dt[, c("name", "batch") := tstrsplit(sample, "/", 2)]
  dt[,sample:=NULL]
  dt$name = paste0(str_trim(dt$name),extension)
  
  location_translation = fread(trans_file, header = TRUE)
  location_translation$Barcode = paste0("barcode",
                                        str_trim(location_translation$Barcode),
                                        extension)
  
  location_translation$`concentration` = as.numeric(location_translation$`concentration`)
  
  dt = merge(dt,location_translation,
             by.x=c('name'),
             by.y=c('Barcode'))
  
  dt$Place = str_trim(dt$Place)
  dt$Date = str_trim(dt$Date)
  dt$Season = str_trim(dt$Season)
  
  
  return(dt)
  
  
}


mybreaks <- function(x) {
  x <- breaks_pretty()(x)
  x[x==0] <- NA
  x
}

global_limits <- c(1, 10^5)

coverage_plotting = function(dt){
  
  dt$Place <- str_trim(gsub("^[A-Z]+", "", dt$Place))
  dt$Place_f = factor(dt$Place, levels=c('Chur','Werdhoelzli','Lugano','Aire'))
  
  
  p <- dt %>%
    ggplot(aes(x=Position,y=Coverage),group = Date)+
    geom_line(stat = 'identity',aes(colour = Date), alpha=0.8,linewidth = 1)+
    scale_y_log10("Read Depth",
                  breaks = trans_breaks("log10", function(x) 10^x[x != 0], n=3),
                  labels = trans_format("log10", math_format(10^.x))
                  
    )+
    scale_x_continuous("Position", 
                       breaks = mybreaks,
                       n.breaks = 3,
                       guide = guide_axis(angle = 30))+
    theme(legend.position="none",
          # title = element_text(size = 16),
          axis.title.y = element_text(margin = margin(r=10), size = 32),
          axis.text.y = element_text(size = 30),
          axis.text.x = element_text(size = 30),
          axis.title.x = element_text(margin = margin(t=10), size = 32),
          strip.text.x = element_text(size = 32),
          strip.text.y = element_text(size = 27.5),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          linewidth = 0.5),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgrey"), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "lightgrey", 
                                      fill=NA, 
                                      linewidth=1))+
    facet_grid(vars(Place_f),vars(Season)) 
  
  return(list(plot = p, plot_data = dt))
  
  
}



amplicon_plotting = function(dt){
  
  p = dt %>% 
    ggplot(aes(x=amplicon,y=read_count,fill= segment)) + 
    geom_boxplot(alpha=0.45,width=0.5,outlier.size = 4)+
    geom_jitter(size=2.7, alpha=0.5,aes(colour = segment),
                position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5))+
    facet_grid(cols = vars(Season))+
    scale_y_log10("Read Count",
                  breaks = trans_breaks("log10", function(x) 10^x[x != 0], n=3),
                  labels = trans_format("log10", math_format(10^.x)),
                  limits = global_limits
    )+ 
    scale_x_discrete(labels= c('1','2','3','4','5','6'),name = "Amplicon")+
    theme(legend.position="bottom",
          legend.title = element_blank(),
          legend.text = element_text(size = 30),
          legend.key.size = unit(2, "lines"),
          axis.title.x = element_text(size = 32, margin = margin(t=15)),
          axis.title.y = element_text(size = 32, margin = margin(r=15)),
          axis.text.x = element_text(size = 30),
          axis.text.y = element_text(size = 30),
          plot.title = element_text(hjust = 0.5),
          strip.text = element_text(size = 32),
          panel.background = element_rect(fill = "white",
                                          colour = "white",
                                          size = 0.5),
          panel.grid.major = element_line(size = 0.25, linetype = 'solid',
                                          colour = "lightgrey"), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(colour = "lightgrey", 
                                      fill=NA, 
                                      linewidth=1)) 
  
  return(p)
  
}







