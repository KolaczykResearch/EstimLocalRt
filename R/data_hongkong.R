library(readr)
library(dplyr)
library(data.table)
library(tidyr)

# set a directory for input
hdir <- "../../Data/HongKong/"

# set a destination for output
sdir <- "../../Data/HongKong/"

# data is available at https://doi.org/10.1038/s41591-020-1092-01
name_csv <- "case_data.csv"  
dir_csv <- paste0(hdir,name_csv)

# daily local and imported infection counts
infection_counts_obs  <-  fread(dir_csv, header = T) %>% select(epi.date,cluster.category) %>%
  mutate(dates=as.Date(epi.date,format="%d/%m/%Y")) %>%
  mutate(incid=ifelse(cluster.category %in% c("Sporadic imported" ,"Cluster of imported cases"),
                      "imported","local")) %>%
  group_by(dates,incid) %>% summarise(count=n()) %>% ungroup %>%
  spread(incid, count)  %>% complete(dates=seq.Date(min(dates), max(dates), by="day"),
                                     fill = list(imported=0,local=0)) %>%
  mutate(dates=format(dates,"20%Y-%m-%d"))

write.csv(infection_counts_obs,file = paste0(sdir,"infection_counts_obs.csv"),row.names = FALSE)   
write.table(c(infection_counts_obs%>%select(imported) %>% unlist/9*10,0),file = paste0(sdir,"imported_counts_adjusted.csv"),row.names = FALSE,col.names = FALSE)   

