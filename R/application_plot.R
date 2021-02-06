library(readr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(data.table)
library(grid)
library(gridExtra)
library(lemon)
 
plot_all <- list() 

####################### Figure 4(a) and (b) #######################
hdir <- "../../Data/HongKong/"
name_csv <- "infection_counts_obs.csv"
dir_csv <- paste0(hdir,name_csv)
infection_counts_obs <- fread(dir_csv, header = T) %>%gather(variables,means,`local`:`imported`,factor_key = TRUE)  %>% 
               mutate(dates=as.Date(dates)) 
plot_date_hongkong <- c(as.Date("2020-01-18"), as.Date("2020-04-25"))
Rt_date_hongkong <- seq(as.Date("2020-01-20"), as.Date("2020-04-25"),by="day")
plot_all[[1]] <- infection_counts_obs %>% mutate(variables=case_when(variables=="local" ~ "Local cases",
                                                                     variables=="imported" ~ "Imported cases"))%>%
  ggplot( aes(as.Date(dates), y = means,fill=variables  ))+ 
  geom_bar(stat="identity", color=NA, position="stack")  +
  xlab("") +ylab("") +ggtitle("Hong Kong: daily infections")+ ylab("")+  ylim(0,50)+
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =plot_date_hongkong) +
  theme_minimal(base_size = 12) + 
  labs(tag = "(a)") +
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())  

hdir <- "../../Data/Australia/"
name_csv <- "infection_counts_obs.csv"
dir_csv <- paste0(hdir,name_csv)
infection_counts_obs <- fread(dir_csv, header = T) %>%gather(variables,means,`local`:`imported`,factor_key = TRUE)  %>% 
      mutate(dates=as.Date(dates)) 
plot_date_australia <- c(as.Date("2020-01-25"), as.Date("2020-04-14"))
Rt_date_australia <- seq(as.Date("2020-03-09"), as.Date("2020-04-14"),by="day")

plot_all[[2]]<- infection_counts_obs %>% mutate(variables=case_when(variables=="local" ~ "Local cases",
                                                                      variables=="imported" ~ "Imported cases"))%>%
  ggplot( aes(as.Date(dates), y = means,fill=variables  ))+ 
  geom_bar(stat="identity", color=NA, position="stack")  +
  xlab("") +ylab("") +ggtitle("Victoria: daily infections")+ ylab("")+  
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =plot_date_australia) +
  theme_minimal(base_size = 12) +  
  labs(tag = "(b)") +
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank())  


####################### Figure 4(c) and (d) #######################
hdir <- "../../Results/Rt_Estimation/"

int_name <- "hongkong"
name_csv <- paste0("application_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
res_Rt_hongkong <-fread(dir_csv, header = T) 

plot_all[[3]]<-res_Rt_hongkong  %>% mutate(Type = recode_factor(Type, `No misidentification` = "alpha[0]==0~`,`~alpha[1]==0 ",
                                                                `Error 1` = "alpha[0]~`~ Beta(2,18),`~alpha[1]~`~ Beta(4,8)` ",
                                                                `Error 2` = "alpha[0]~`~ Beta(4,8),`~alpha[1]~`~Beta(2,18)` " )) %>%
  filter(dates %in% as.character(Rt_date_hongkong))%>%
  ggplot(aes(as.Date(dates), y = means,color=Type) ) +
  geom_ribbon( aes(ymin = Low, ymax = High, fill = Type, color = NULL), alpha = .15) +
  geom_line( aes(y = means), size = 1)+ 
  xlab("") + ylab("") +   ylim(0,6)+ ggtitle("Hong Kong: local time-vary reproduction numbers")+
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =  plot_date_hongkong) +
  theme_minimal(base_size = 12) +  
  labs(tag = "(c)") +
  geom_hline(yintercept=1, linetype="dashed")+
  scale_colour_discrete(labels = parse_format()) +
  scale_fill_discrete(labels = parse_format())  +
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank(),
        legend.text.align = 0)

int_name <- "australia"
name_csv <- name_csv <- paste0("application_",int_name,".csv") 
dir_csv <- paste0(hdir,name_csv)
res_Rt_australia <-fread(dir_csv, header = T) 

plot_all[[4]]<- res_Rt_australia  %>% mutate(Type = recode_factor(Type, `No misidentification` = "alpha[0]==0~`,`~alpha[1]==0 ",
    `Error 1` = "alpha[0]~`~ Beta(2,18),`~alpha[1]~`~ Beta(4,8)` ",
    `Error 2` = "alpha[0]~`~ Beta(4,8),`~alpha[1]~`~Beta(2,18)` " ))%>%
  filter(dates %in% as.character(Rt_date_australia))%>%
  ggplot(aes(as.Date(dates), y = means,color=Type) ) +
  geom_ribbon( aes(ymin = Low, ymax = High, fill = Type, color = NULL), alpha = .15  ) +
  geom_line( aes(y = means), size = 1)+ 
  xlab("") + ylab("") +   ylim(0,6)+ ggtitle("Victoria: local time-vary reproduction numbers")+
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =  plot_date_australia) +
  theme_minimal(base_size = 12) + 
  labs(tag = "(d)") +
  geom_hline(yintercept=1, linetype="dashed")+
  scale_colour_discrete(labels = parse_format()) +
  scale_fill_discrete(labels = parse_format())  +
theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank(),
      legend.text.align = 0)
 
