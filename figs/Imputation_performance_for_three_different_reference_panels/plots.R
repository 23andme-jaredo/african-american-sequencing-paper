library(ggplot2)
w <- 4
bakeoff <- read.csv('1kgp-hrc-afam-bakeoff.csv',as.is=TRUE)
bakeoff$variant <- factor(toupper(bakeoff$variant),levels=c('SNP','INDEL'))
bakeoff$panel <- factor(gsub('DV-0.10.0','AFAM-DV-GLx',gsub('1000G','1KGP',bakeoff$panel)),levels=c("AFAM-DV-GLx","1KGP","HRC"))

m <- 1
min.af <- 1e-4
plt <- ggplot(bakeoff,aes(x=af,y=r2,col=panel,lty=variant)) +geom_line() +theme_bw()
plt <- plt+scale_color_brewer(palette='Set1')+xlab("gnomAD AFR frequency")
plt <- plt+scale_x_log10(limits=c(min.af,1),labels=function(x) format(x, big.mark = ",", scientific = FALSE))
plt <- plt+theme(legend.title=element_blank(),legend.position=c(.01,.995),legend.justification=c(0,1),legend.margin = margin(0,.1,.1,0))

plt <- plt+scale_y_continuous(breaks=seq(0,1,.1),limits=c(0,1),name=bquote('Aggregate '~R^2))
ggsave('three-panel-bakeoff.png',plot=plt,width=w,height=w)

for(p in unique(bakeoff$panel)) message(p,' = ',with(subset(bakeoff,variant=='SNP' & panel==p),approx(af,r2,0.005))$y)

plotme <- read.csv('gnomad-sensitivity.csv',as.is=T)
plotme$variant <- factor(toupper(plotme$variant),levels=c('SNP','INDEL'))

rename <- c('HRC','1KGP','AFAM-DV-GLx')
names(rename) <- c('hrc','ogp','afam')
plotme$panel <- factor(rename[plotme$name],levels=c("AFAM-DV-GLx","1KGP","HRC"))

plt <- ggplot(plotme,aes(x=gnomad,y=value,col=panel,lty=variant)) + geom_line()+theme_bw()+scale_color_brewer(palette='Set1') +scale_x_log10(limits=c(min.af,1),labels=function(x) format(x, big.mark = ",", scientific = FALSE))+scale_y_continuous(breaks=seq(0,1,.1))+theme(legend.position="none",)+xlab("gnomAD AFR frequency")+ylab("Proportion of gnomAD variants present in panel")
ggsave('gnomad-sensivity.png',plot=plt,width=w,height=w)

for(p in rename) message(p,' = ',with(subset(plotme,variant=='SNP' & panel==p),round(100*approx(gnomad,value,0.005)$y,1)))


