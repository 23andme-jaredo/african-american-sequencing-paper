library(ggplot2)
#library(gridExtra)
#library(kable)

dat <- read.csv('single-sample-happy-evaluation.csv',as.is=TRUE)
dat1 <- subset(dat,Type=='SNP')

plt <- ggplot(dat1,aes(x=coverage,y=METRIC.F1_Score,col=caller))+ theme_bw()+ geom_point(alpha=.9)+geom_abline(intercept=0,slope=1) + theme(legend.position=c(.95,0.05),legend.justification=c(1,0),legend.title=element_blank())  +scale_color_brewer(palette='Set1') +xlab("Coverage") + ylab("SNP F1")
plt

ggsave('f1-versus-depth.png',plt,width=5,height=5)
