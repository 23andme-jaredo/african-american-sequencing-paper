library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
w <- 4

cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")[c(7,6,4,8,2)]

make.plot <- function(bakeoff) {
    bakeoff$panel <- factor(bakeoff$panel,levels=c("AFAM-DV-GLx","1KGP","HRC",'TOPMed','CAAPA'))
    bakeoff$variant <- factor(bakeoff$variant,levels=c('SNP','INDEL'))
    m <- 1
    min.af <- 1e-5
    plt <- ggplot(bakeoff,aes(x=af,y=r2,col=panel,lty=variant)) +geom_line() +theme_bw()
    plt <- plt+xlab("gnomAD AFR frequency")
    plt <- plt+scale_colour_manual(values=cbPalette)
    plt <- plt +theme(legend.position='none')
    plt <- plt+scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels=trans_format("log10", math_format(10^.x)))
    plt <- plt+scale_y_continuous(breaks=seq(0,1,.1),limits=c(0,1),name=bquote('Aggregate '~R^2))
    plt
}

for(m in c('beagle','minimac4')){
    bakeoff <- subset(read.csv('5-way-bakeoff.csv',as.is=TRUE),method==m)
    plt <- make.plot(bakeoff)+theme(legend.title=element_blank(),legend.position=c(.995,0.005),legend.justification=c(1,0),legend.margin = margin(0,.1,.1,0))
    ggsave(paste0(m,'-five-panel-bakeoff.png'),plot=plt,width=w,height=w,dpi=600)

    bakeoff <- subset(read.csv('5-way-bakeoff-hard.csv',as.is=TRUE),method==m)
    ggsave(paste0(m,'-five-panel-bakeoff-hard.png'),plot=make.plot(bakeoff),width=w,height=w,dpi=600)
}

bakeoff <- subset(read.csv('5-way-bakeoff-hard.csv',as.is=TRUE),method=='minimac4')

tab <- data.frame(panel=unique(bakeoff$panel,r2=NA))
for(i in 1:nrow(tab)) tab$r2[i] <- with(subset(bakeoff,variant=='SNP' & panel==tab$panel[i]),round(approx(af,r2,0.005)$y,2))
tab

fname <- 'gtex-sensitivity.csv'
rename <-        c('AFAM-DV-GLx','1KGP','HRC','TopMED','CAAPA')
names(rename) <- c('afam',       'ogp' ,'hrc','topmed','caapa')
plotme <- pivot_longer(read.csv(fname),cols=names(rename),names_to='panel')
plotme$panel <- factor(rename[plotme$panel],levels=rename)
plotme$bin <- with(plotme,factor(bin,levels=unique(bin)))
plotme$pct <- round(100*plotme$value)
plotme$variant <- factor(plotme$variant,levels=c('SNP','INDEL'))
plotme <- subset(plotme,!(panel%in%c('HRC','CAAPA') & variant=='INDEL'))

plt <- ggplot(plotme,aes(x=as.integer(bin),y=100*value,col=panel,lty=variant))
plt <- plt+geom_point()+theme_bw()+scale_colour_manual(values=cbPalette)
plt <- plt+geom_line()+xlab("GTEx allele count")+ylab("% present in panel")
plt <- plt+scale_y_continuous(breaks=seq(0,100,10))+theme(legend.position='none')
plt <- plt+scale_x_continuous(labels=levels(plotme$bin))
                                        #plt <- make.plot2('gtex-sensitivity.csv')
ggsave('gtex-sensivity.png',plot=plt,width=w,height=w,dpi=600)

#system('convert +append minimac4-five-panel-bakeoff-hard.png minimac4-five-panel-bakeoff.png gtex-sensivity.png Figure4.png')
arrange(pivot_wider(plotme,id_cols=c('panel','variant','bin'),values_from='pct',names_from='panel'),variant,bin)
#for(p in rename) message(p,' = ',with(subset(plotme,variant=='SNP' & panel==p),round(100*approx(gtex.ac,value,0.005)$y,1)))
