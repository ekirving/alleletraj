library(dplyr)
library(stringr)
library(ggplot2)

setwd('/Users/Evan/Dropbox/Code/alleletraj')

dfr <- read.delim("./bed/eud-thin.ld.summary",sep="",header=F,check.names=F,stringsAsFactors=F)
colnames(dfr) <- c("dist","rsq")

dfr$distc <- cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=10000))

dfr1 <- dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))

dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))

dist_limit <- 10^6

ggplot()+
    geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
    geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
    labs(x="Distance (Kb)",y=expression(LD~(r^{2})))+
    scale_x_continuous(breaks=seq(0, dist_limit, dist_limit/10),labels=seq(0, dist_limit, dist_limit/10)/1000, limits = c(0,dist_limit))+
    theme_bw()
