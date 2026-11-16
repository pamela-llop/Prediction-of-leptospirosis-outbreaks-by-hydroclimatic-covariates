
# Plots ------------------------

# library that are needed
library(ggsci)
library(reshape2)
library(ggplot2)
library(cowplot)

# plot of variables -----------------------------

# load tha data: copy the dataset file into CODES folder
# data has four columns: 
# column 1: the response (leptospirosis cases)
# columns 2 to 4: covariates (precipitation, hydrometric level, ONI)

source("~/CODES/data")

t=c(1:144)
df1=data.frame(t,data[,1])
df2=data.frame(t,data[,2])
df3=data.frame(t,data[,3])
df4=data.frame(t,data[,4])
etiquetas=c('2009', '2010', '2011', '2012', '2013', '2014', '2015', '2016', '2017', '2018', '2019','2020')

g1=ggplot(df1,aes(t,data[,1]))+geom_line()+labs(x='',y="Leptospirosis Cases")+
  ggtitle ("City name") +  
  theme_bw() +
  theme (plot.title = element_text(size=25,vjust=1, face="bold", lineheight=1.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_x_continuous(breaks=seq(1,N,12), labels=etiquetas)+ 
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 12))
g2=ggplot(df2,aes(t,data[,2]))+geom_line()+labs(x='',y="Precipitation")+
  theme_bw() +
  scale_x_continuous(breaks=seq(1,N,12), labels=etiquetas)+   
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 12))
g3=ggplot(df3,aes(t,data[,3]))+geom_line()+labs(x='',y="Hydrometric Level")+
  theme_bw() +
  scale_x_continuous(breaks=seq(1,N,12), labels=etiquetas)+
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 12))
g4=ggplot(df4,aes(t,data[,4]))+geom_line()+labs(x= 'Year',y="ONI")+
  theme_bw() +
  scale_x_continuous(breaks=seq(1,dim(data)[1],12), labels=etiquetas)+
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 12))+
  theme(axis.title.x = element_text(face="bold", vjust=1.2, size=15))


plot_grid(g1,g2,g3,g4,ncol = 1)
ggsave("variables.jpeg", width = 13, height = 30, units = "cm") 


# plots of annual curvs of leptospirosis cases ------------------------

x <- t(matrix(data[,1], ncol = 12))[1:12,] 
rownames(x) <- c('2009','2010','2011','2012','2013','2014','2015','2016','2017','2018','2019','2020')
colnames(x) <- 1:12 
y <- melt(x)
colnames(y) <- c("Year", "Month", "Cases")
ggplot(y) + geom_line(aes(x=Month, y=Cases, group=Year), size = 1) + 
  theme_bw() +
  labs(x= 'Month',y="Number of Cases") +
  ggtitle ("City name") +  
  theme (plot.title = element_text(size=25,vjust=1, face="bold", lineheight=1.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.title.x = element_text(face="bold", vjust=1.2, size=18))+ 
  theme(axis.text.x = element_text(vjust = 1, hjust=1, size = 10))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 10))+
  scale_x_continuous(breaks=1:12, labels=1:12)+
  facet_wrap(~Year, nrow = 2)+
  theme(strip.text.x = element_text(face = "bold", size=13))

ggsave("curvs.jpeg", width = 35, height = 10, units = "cm") 


# after running Analysis.R, save the results in results.Data and then run:

load('results.RData')
pred_SFPLR=round(save_list[[1]])
ar_pred=round(save_list[[2]])
arx_pred=round(save_list[[3]])
kpred=save_list[[4]]
kkpred=save_list[[5]]
real=save_list[[6]]


# plots of predictions ----------------------------

x=rbind(real,round(ar_pred),round(arx_pred),
        round(pred_SFPLR),round(kpred),round(kkpred)) 
colnames(x) <- factor(1:12)
rownames(x) <- c("Observed", "ARIMA", "ARIMAX", "SFPLR", "NW1K", "NW2K") 
y <- melt(x)
colnames(y) <- c("Curv", "Month", "Cases")
ggplot(y) + geom_line(aes(x=Month, y=Cases, group=Curv), size = 1) +
  theme_bw() +
  labs(x= 'Month',y="Number of Cases") +
  ggtitle ("City name") + 
  theme (plot.title = element_text(size=25,vjust=1, face="bold", lineheight=1.5))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(face="bold", vjust=1.2, size=18))+
  theme(axis.title.x = element_text(face="bold", vjust=1.2, size=18))+ 
  theme(axis.text.x = element_text(vjust = 1, hjust=1, size = 10))+
  theme(axis.text.y = element_text(vjust = 1, hjust=1, size = 10))+
  scale_x_continuous(breaks=1:12, labels=1:12)+
  facet_wrap(~Curv, nrow = 2)+
  theme(strip.text.x = element_text(face = "bold", size=15))
ggsave("predictions.jpeg", width = 18, height = 15, units = "cm")


