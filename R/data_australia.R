library(readr)
library(dplyr)
library(data.table)
library(readxl)

# set a directory for input
hdir <- "../../Data/Australia/"

# set a destination for output
sdir <- "../../Data/Australia/"

# data is available at https://doi.org/10.1101/2020.05.12.20099929
name_csv <- "41467_2020_18314_MOESM4_ESM.xlsx"
dir_csv <- paste0(hdir,name_csv)
data1  <-  read_excel(dir_csv)  %>% select(VIC_ID,Overseas_travel)%>% filter(!is.na(Overseas_travel))
name_csv <- "VIC.csv"
dir_csv <- paste0(hdir,name_csv)
data2 <-  fread(dir_csv,header = T)  %>% select(VIC_ID,Date_coll) %>% filter(!is.na(Date_coll))

# daily local and imported infection counts
infection_counts_obs  <- merge(data1,data2,by="VIC_ID") %>%
  mutate(dates=as.Date(Date_coll,format="%d/%m/%Y")) %>%
  mutate(incid=ifelse(Overseas_travel == "Yes", "imported","local")) %>%
  group_by(dates,incid) %>% summarise(count=n()) %>% ungroup %>%
  spread(incid, count) %>% complete(dates=seq.Date(min(dates), max(dates), by="day"),
                                    fill = list(imported=0,local=0))%>%
  mutate(dates=format(dates,"%Y-%m-%d"))

write.csv(infection_counts_obs,file = paste0(sdir,"infection_counts_obs.csv"),row.names = FALSE)   
write.table(c(infection_counts_obs%>%select(imported) %>% unlist/9*10,rep(0,7)),file = paste0(sdir,"imported_counts_adjusted.csv"),row.names = FALSE,col.names = FALSE)   

 
