library(scales)
library(ggplot2)
library(data.table)

bakeoff <- read.csv('gatk-versus-dv-panel-evaluation.csv'
                   ,as.is=TRUE
                   ,col.names=c('type','count','AF','r2','panel','truth')
                   ,header=FALSE)

bakeoff$type <- factor(bakeoff$type,levels=c('SNP','INDEL'))
cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(7,6,4)]

plt <- ggplot(bakeoff,aes(x=AF,y=r2,col=panel,shape=panel))
plt <- plt+geom_line()+theme_bw()+geom_point()
plt <- plt+xlab("Alternate allele frequency")+scale_colour_manual(values=cbPalette)
plt <- plt+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels=trans_format("log10", math_format(10^.x)))
plt <- plt+scale_y_continuous(breaks=seq(0,1,.1),limits=c(0.2,1),name=bquote('Aggregate '~R^2))
plt <- plt+facet_grid(cols=vars(truth),rows=vars(type))+theme(legend.position=c(0.4,0.1),legend.title=element_blank())
plt <- plt+scale_shape_manual(values = c(16,18))

lbls <- expand.grid(type=c('SNP','INDEL'),truth=c('DV truth','GATK truth'))
lbls$AF <- .0001
lbls$r2 <- 1
lbls$txt <- c('a','c','b','d')
lbls$panel <- NA
plt <- plt+geom_text(aes(label=txt),lbls,col='black',size=5)
plt

ggsave('Figure3.png',width=8,height=8,dpi=600)
