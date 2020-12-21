library(ggplot2)
library(data.table)
bakeoff <- read.csv('../../data/gatk-versus-dv-panel-evaluation.csv',as.is=TRUE)
bakeoff <- cbind(bakeoff[,1:4],(setNames(data.frame(tstrsplit(bakeoff$panel,'with')),c('panel','truth'))))
bakeoff$type <- factor(bakeoff$type,levels=c('SNP','INDEL'))

plt <- ggplot(bakeoff,aes(x=AF,y=r2,col=panel))+geom_line()+scale_x_log10()+theme_bw()+geom_point()+scale_color_brewer(palette='Set1')+xlab("Alternate allele frequency")
plt+scale_y_continuous(breaks=seq(0,1,.1),limits=c(0.2,1),name=bquote('Aggregate '~R^2))+facet_grid(rows=vars(truth),cols=vars(type))+theme(legend.position=c(0.4,0.1))
ggsave('gatk-versus-dv.png',width=8,height=8)
