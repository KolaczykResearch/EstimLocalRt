library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(data.table)
library(grid)
library(gridExtra)
library(lemon)

# set a directory for input
hdir <- "../../Results/Simulation_Epidemics/"


# mean of daily local and imported diagnosed counts in 1,000 simulation trials forepidemics in Hong Kong
int_name <- "hongkong"
date_seq <- seq(as.Date("2020-01-18"), as.Date("2020-04-25"), by = "day")
name_csv <- paste0("diagnosed_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
n_sims <- 1000
sim_results_mean_diagnosed1 <-  fread(dir_csv, header = T) %>% select(sim_num,dates,exogenous)%>%
  mutate(exogenous=case_when(exogenous ==1 ~ "imported",exogenous== 0 ~ "local"))    %>%
  group_by(sim_num,dates,exogenous) %>% summarise(count=n())   %>%
  ungroup %>% complete(sim_num= 0:(n_sims-1),dates= as.character(date_seq), exogenous,fill = list(count = 0))%>%
  mutate(dates=as.Date(dates))%>%group_by(dates,exogenous) %>% summarise(means=mean(count))%>% 
  spread(exogenous,  means) %>%gather(variables,means,`local`:`imported`,factor_key = TRUE) %>% mutate(country="Hong Kong")

# mean of daily local and imported diagnosed counts in 1,000 simulation trials forepidemics in Victoria, Australia
int_name <- "australia"
date_seq <- seq(as.Date("2020-01-25"), as.Date("2020-04-14"), by = "day")
name_csv <- paste0("diagnosed_counts_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
sim_results_mean_diagnosed2 <-  fread(dir_csv, header = T) %>% select(sim_num,dates,exogenous)%>%
  mutate(exogenous=case_when(exogenous ==1 ~ "imported",exogenous== 0 ~ "local"))    %>%
  group_by(sim_num,dates,exogenous) %>% summarise(count=n())   %>%
  ungroup %>% complete(sim_num= 0:(n_sims-1),dates= as.character(date_seq), exogenous,fill = list(count = 0))%>%
  mutate(dates=as.Date(dates))%>%group_by(dates,exogenous) %>% summarise(means=mean(count))%>% 
  spread(exogenous,  means) %>%gather(variables,means,`local`:`imported`,factor_key = TRUE) %>% mutate(country="Victoria") %>% 
  filter(dates %in% date_seq)

 
plot_diagnosed <- list()
plot_diagnosed[[1]]  <- sim_results_mean_diagnosed1%>% mutate(variables=case_when(variables=="local" ~ "Local cases",
   variables=="imported" ~ "Imported cases")) %>% ggplot( aes(as.Date(dates), y = means,fill=variables  ))+ 
  geom_bar(stat="identity", color=NA, position="stack")  +
  xlab("") +ylab("") +ggtitle("Daily diagnosed counts in Hong Kong")+ ylab("")+  
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d") + ylim(0,50)+
  theme_minimal(base_size = 12) +  
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())  

plot_diagnosed[[2]]  <- sim_results_mean_diagnosed2%>% mutate(variables=case_when(variables=="local" ~ "Local cases",
                                                                                  variables=="imported" ~ "Imported cases")) %>% ggplot( aes(as.Date(dates), y = means,fill=variables  ))+ 
  geom_bar(stat="identity", color=NA, position="stack")  +
  xlab("") +ylab("") +ggtitle("Daily diagnosed counts in Victoria")+ ylab("")+  
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d") + ylim(0,50)+
  theme_minimal(base_size = 12) +  
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())  

 

