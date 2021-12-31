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
hdir <- "../Results/Rt_Estimation/"
 

# plots for estimations of local time-varying reproduction numbers in simulated epidemics for Hong Kong under different error misidentification rates
int_name <- "hongkong"
n_alpha <- 6
name_csv <- paste0("simulation_",int_name,"_", 1:n_alpha,".csv") 
dir_csv <- paste0(hdir,name_csv)
res_Rt_hongkong <- lapply(1:n_alpha, function(i) fread(dir_csv[i], header = T) )   
 
ymax <- c(rep(4,3),rep(5,3))
plot_date_hongkong <- c(as.Date("2020-01-18"), as.Date("2020-04-25"))
Rt_date_hongkong <- seq(as.Date("2020-01-25"), as.Date("2020-04-25"),by="day")

mean_alpha_0 <- sapply(1:n_alpha, function(i) res_Rt_hongkong[[i]]$zeta_0_true[1]/(res_Rt_hongkong[[i]]$xi_0_true[1]+res_Rt_hongkong[[i]]$zeta_0_true[1]) )
mean_alpha_1 <- sapply(1:n_alpha, function(i) res_Rt_hongkong[[i]]$zeta_1_true[1]/(res_Rt_hongkong[[i]]$xi_1_true[1]+res_Rt_hongkong[[i]]$zeta_1_true[1]) )

plot_Rt_hongkong <- lapply(1:n_alpha, function(i)  res_Rt_hongkong[[i]] %>% filter(dates %in% as.character(Rt_date_hongkong))%>% 
                             filter(Type != "Misidentified_EpiEstim(no sliding)")%>%
ggplot(aes(as.Date(dates), y = means,color=Type) ) +
  geom_ribbon( aes(ymin = Low, ymax = High, fill = Type, color = NULL), alpha = .15) +
  geom_line( aes(y = means), size = 1)+ 
  xlab("") + ylab("Local reproduction numbers") +   ylim(0,ymax[i])+
  ggtitle(bquote("Hong Kong :"~alpha[0]~"~"~list(.(mean_alpha_0[i]))~","~alpha[1]~"~"~list(.(mean_alpha_1[i]))))+
  scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =  plot_date_hongkong) +
  theme_minimal(base_size = 16) +  
  geom_hline(yintercept=1, linetype="dashed")+
  theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=0.5),
        panel.background = element_blank()) )


 
# plots for estimations of local time-varying reproduction numbers in simulated epidemics for Victoria, Australia under different error misidentification rates
int_name <- "australia"
name_csv <- paste0("simulation_",int_name,"_", 1:n_alpha+n_alpha,".csv") 
dir_csv <- paste0(hdir,name_csv)
res_Rt_australia <- lapply(1:n_alpha, function(i) fread(dir_csv[i], header = T) )   
plot_date_australia <- c(as.Date("2020-01-25"), as.Date("2020-04-14"))
Rt_date_australia <- seq(as.Date("2020-03-09"), as.Date("2020-04-14"),by="day")

plot_Rt_australia <- lapply(1:n_alpha, function(i)  res_Rt_australia[[i]] %>% filter(dates %in% as.character(Rt_date_australia))%>%
                             filter(Type != "Misidentified_EpiEstim(no sliding)")%>%
                             ggplot(aes(as.Date(dates), y = means,color=Type) ) +
                             geom_ribbon( aes(ymin = Low, ymax = High, fill = Type, color = NULL), alpha = .15) +
                             geom_line( aes(y = means), size = 1)+ 
                             xlab("") + ylab("Local reproduction numbers") +   ylim(0,ymax[i])+
                             ggtitle(bquote("Victoria :"~alpha[0]~"~"~list(.(mean_alpha_0[i]))~","~alpha[1]~"~"~list(.(mean_alpha_1[i]))))+
                             scale_x_date(date_breaks = "15 day", date_labels =  "%b-%d",limits =plot_date_australia  ) +
                             theme_minimal(base_size = 16) +  
                             geom_hline(yintercept=1, linetype="dashed")+
                             theme(legend.title = element_blank(),legend.position = c(0.15, 0.85),
                                   plot.title = element_text(hjust = 0.5),
                                   panel.grid.minor = element_blank(),
                                   panel.border = element_rect(fill=NA,color="black", size=0.5),
                                   panel.background = element_blank()) )

sdir <- "../Results/Figure/"
pdf(file=paste0(sdir,"Fig3.pdf"), width =12, height =12,onefile = F)
grid_arrange_shared_legend( plot_Rt_hongkong [[1]],plot_Rt_australia[[1]] ,plot_Rt_hongkong [[2]],plot_Rt_australia[[2]] ,
                            plot_Rt_hongkong [[3]], plot_Rt_australia [[3]], ncol = 2, nrow = 3, position='bottom')
dev.off()

pdf(file=paste0(sdir,"Fig4.pdf"), width =12, height =12,onefile = F)
grid_arrange_shared_legend( plot_Rt_hongkong [[4]],plot_Rt_australia[[4]] ,plot_Rt_hongkong [[5]],plot_Rt_australia[[5]] ,
                            plot_Rt_hongkong [[6]], plot_Rt_australia [[6]], ncol = 2, nrow = 3, position='bottom')
dev.off()


 