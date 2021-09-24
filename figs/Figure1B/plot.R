library(ggplot2)
library(Hmisc)

df <- read.csv('ancestry_proportions.csv',as.is=TRUE,header=TRUE)
exclude <- c("European", "Native American/East Asian")
lvls <- with(df,c(setdiff(pretty.leaf[order(pct.all)],exclude),exclude))
df$pretty.leaf <- with(df,factor(pretty.leaf,levels=lvls))
lvls <- with(df,round(tapply(pct.all,pretty.node,sum),1))
df$node.pct.all <- with(df,lvls[pretty.node])
lvls <-  with(df,c(setdiff(unique(pretty.node[order(node.pct.all,decreasing=TRUE)]),exclude),exclude))
df$pretty.node2 <- with(df,paste0(pretty.node,": ",node.pct.all,'%'))
lvls <- with(subset(df,!duplicated(pretty.node2)),c(setdiff(pretty.node2[order(node.pct.all)],pretty.node2[pretty.node%in%exclude]),pretty.node2[pretty.node%in%exclude]))
df$pretty.node2 <- with(df,factor(pretty.node2,levels=lvls))

m <- 40
m2 <- 5

myplot <- function(my.y) {
    cbPalette <- c("#D55E00", "#56B4E9",  "#999999","#E69F00", "#009E73","#F0E442", "#0072B2", "#CC79A7","#377EB8","#4DAF4A")
    cbPalette <- c("#D55E00", "#CC79A7", "#999999","#F0E442","#E69F00","#4DAF4A","#377EB8")
    plt <- ggplot(df,aes_string(x='pretty.leaf',y=my.y,fill='pretty.node2'))+geom_bar(stat='identity')+theme_bw()
    plt <- plt+xlab('')+theme(axis.text.x = element_text(angle=45,hjust=1))
    plt <- plt+ylab("% Ancestry contribution")
    plt <- plt+scale_fill_manual(values=cbPalette)
    plt#+theme(legend.position='none')
}

plt <- myplot('pct.all')
plt

#plt <- plt+theme(legend.position=c(0.02,.98),legend.justification=c(0,1),legend.title=element_blank(),plot.margin=margin(t=m2,r=m2,b=m2,l=m, unit = "pt"))
plt <- plt+theme(legend.title=element_blank(),plot.margin=margin(t=m2,r=m2,b=m2,l=m, unit = "pt"))
plt+annotate("text", -Inf, Inf, label = "Top-left", hjust = 0, vjust = 1)
plt

ggsave("Figure_1B.png",plt,width=8,height=5,dpi=600)
#for(i in 1:6) ggsave(paste0('cluster',i,"-barplot.png"),myplot(paste0('pct.cluster',i)),width=5,height=5)

tab <- df[order(df$node.pct,decreasing=TRUE),c('pretty.node','pretty.leaf',paste0('pct.cluster',1:6))]
write.csv(tab,row.names=FALSE)
