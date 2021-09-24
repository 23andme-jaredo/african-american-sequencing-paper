library(ggplot2)
library(tidyr)

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(7,4,2)]

dat <- read.table('Coverage_Titration_Errors.tsv',as.is=T,header=T)
dat <- pivot_longer(dat,cols=!Coverage,names_to='Caller',values_to='Errors')
dat$Caller <- sub('\\.','-',dat$Caller)
plt <- ggplot(subset(dat,Caller!='Freebayes'),aes(x=Coverage,col=Caller,shape=Caller,y=Errors/1000))
plt <- plt+geom_line() +theme_bw()+theme(legend.position=c(.95,.95),legend.justification=c(1,1),legend.title=element_blank())
plt <- plt+ylab("Total errors (1000s)")+geom_point()
plt <- plt+scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))
plt <- plt +scale_colour_manual(values=cbPalette)
plt <- plt#+annotate("text", -Inf, Inf, label = " b", hjust = 0, vjust = 1)
w <- 5
ggsave('Figure_2A.png',plot=plt,width=w,height=w,dpi=600)
