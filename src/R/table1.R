library(tidyr)
library(dplyr)
library(data.table)

dat <- fread('../../figs/f1-versus-depth/single-sample-happy-evaluation.csv')

build.table <- function(dataset) {
    table1 <- summarise(group_by(dataset,Type,caller),
                        N=n(),
                        Recall=sum(TRUTH.TP)/sum(TRUTH.TOTAL),
                        Precision=sum(TRUTH.TP)/(sum(TRUTH.TP)+sum(QUERY.FP)),
                        TP=sum(TRUTH.TP),
                        FN=sum(TRUTH.FN),
                        FP=sum(QUERY.FP),
                        FP.gt=sum(FP.gt),
                        FP.al=sum(FP.al))

    table1$F1 <- with(table1,2*Recall*Precision/(Recall+Precision))
    cn <- c('caller','Type','N','F1','Recall','Precision','TP','FN','FP','FP.gt','FP.al')
    table1 <- table1[,cn]
    print(table1,n=100)
    table1
}

write.csv(build.table(dat),row.names=F)
