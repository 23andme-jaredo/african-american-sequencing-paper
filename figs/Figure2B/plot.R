library(ggplot2)

dat <- read.csv('single-sample-happy-evaluation.csv',as.is=TRUE)
dat$caller <- factor(dat$caller,levels=c("DV-0.10.0","GATK-4.1.0.0","GATK-3.5.0"))

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(7,4,6)]

plt <- ggplot(dat,aes(x=coverage,y=METRIC.F1_Score,col=caller,shape=caller))
plt <- plt+theme_bw()+geom_point(alpha=.9)+geom_abline(intercept=0,slope=1)
plt <- plt+theme(legend.position=c(.95,0.05),legend.justification=c(1,0),legend.title=element_blank())
plt <- plt+xlab("Coverage")+ylab("SNP F1")+scale_colour_manual(values=cbPalette)
plt <- plt+scale_shape_manual(values = c(16,17,18))

ggsave('Figure_2B.png',plt,width=5,height=5,dpi=600)
