library(ggplot2)
library(tidyr)
dat <- read.table('Coverage_Titration_Errors.tsv',as.is=T,header=T)
dat <- pivot_longer(dat,cols=!Coverage,names_to='Caller',values_to='Errors')
plt <- ggplot(subset(dat,Caller!='Freebayes'),aes(x=Coverage,col=Caller,y=Errors/1000)) + geom_line() + theme_bw()+theme(legend.position=c(.95,.95),legend.justification=c(1,1))+ylab("Total errors (1000s)")+geom_point()+ scale_y_continuous(labels=function(x) format(x, big.mark = ",", scientific = FALSE))+scale_color_manual(values=c("#E41A1C","#377EB8","#FF7F00"))
plt
w <- 5
ggsave('num_errors_against_coverage.png',plot=plt,width=w,height=w)

# +scale_color_brewer(palette='Set1')
