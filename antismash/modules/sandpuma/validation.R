setwd("~/git/antismash5/antismash/antismash/modules/sandpuma/")
library(dplyr)
pplacer_scored <- read.table("jk_scored.tsv", sep="\t", header=T) %>% mutate(source='sp2')
original_scored <- read.table("original_jk.tsv", sep="\t", header=T) %>% mutate(source='sp1')
original_scored$query <- gsub("\\(", "-", original_scored$query)
original_scored$query <- gsub("\\)", "-", original_scored$query)
original_scored$query <- gsub("\\.", "-", original_scored$query)



predicat <- original_scored %>% select(shuffle, jackknife, query, pid, spec, called_spec, method, call_made, call, methshuf, source) %>% filter(method=='prediCAT') %>% mutate(snn=NA) %>% rbind(pplacer_scored %>% filter(method=='prediCAT'))
pc.sum.acc <- predicat %>% filter(called_spec!='no_confident_result') %>% group_by(source, method, methshuf, call) %>% summarize(count=n())
predicatf <- original_scored %>% select(shuffle, jackknife, query, pid, spec, called_spec, method, call_made, call, methshuf, source) %>% filter(method=='forced_prediCAT') %>% mutate(snn=NA) %>% rbind(pplacer_scored %>% filter(snn>=0.5 & method=='forced_prediCAT'))
pcf.sum.acc <- predicatf %>% filter(called_spec!='no_confident_result') %>% group_by(source, method, methshuf, call) %>% summarize(count=n())
acc.w <- rbind(pcf.sum.acc, pc.sum.acc) %>% dcast(method+source+methshuf ~ call) %>% mutate(pct=100*(Y/(Y+N)))
acc.w$mdisp <- 'Monophyly'
acc.w[acc.w$method=='forced_prediCAT',]$mdisp <- 'SNN >= 0.5'
acc.w$sdisp <- 'SANDPUMA1 (FastTree)'
acc.w[acc.w$source=='sp2',]$sdisp <- 'SANDPUMA2 (pplacer)'
ggplot(acc.w, aes(x=mdisp, y=pct, color=sdisp))+
  geom_boxplot() + theme_classic() +
  theme(
    legend.position = 'top',
    legend.title = element_blank()
  ) +
  scale_color_brewer(palette='Dark2') +
  ylab("Percent Correct\n(90/10 Cross-validation; 10x)") + xlab("prediCAT Method")


pc.sum.cov <- predicat %>% group_by(source, method, methshuf, call_made) %>% summarize(count=n())
pcf.sum.cov <- predicatf %>% group_by(source, method, methshuf, call_made) %>% summarize(count=n())
cov.w <- rbind(pcf.sum.cov, pc.sum.cov) %>% dcast(method+source+methshuf ~ call_made) 
cov.w[is.na(cov.w)] <- 0
cov.w <- cov.w %>% mutate(pct=100*(Y/(Y+N)))
cov.w$mdisp <- 'Monophyly'
cov.w[acc.w$method=='forced_prediCAT',]$mdisp <- 'SNN >= 0.5'
cov.w$sdisp <- 'SANDPUMA1 (FastTree)'
cov.w[acc.w$source=='sp2',]$sdisp <- 'SANDPUMA2 (pplacer)'
ggplot(cov.w, aes(x=mdisp, y=pct, color=sdisp))+
  geom_boxplot() + theme_classic() +
  theme(
    legend.position = 'top',
    legend.title = element_blank()
  ) +
  scale_color_brewer(palette='Dark2') +
  ylab("Percent Covered\n(90/10 Cross-validation; 10x)") + xlab("prediCAT Method")

snn <- pplacer_scored %>% filter(method=='forced_prediCAT' & called_spec!='no_confident_result')
snn$bin <- as.numeric(as.character(cut(snn$snn, breaks = seq(0,2.5,.5), labels=seq(0,2,.5))))
snn[snn$snn >= 2.5, ]$bin <- 2.5
snn[is.na(snn)] <- 0
snn$callb <- 1
snn[snn$call=='N',]$callb <- 0
bysnn <- snn %>% group_by(bin) %>% summarize(N=length(callb), mean=mean(callb), sd=sd(callb), se=sd/sqrt(N))
bysnn$se[is.na(bysnn$se)] <- 0
bysnn$ymax <- bysnn$mean + bysnn$se
bysnn$ymin <- bysnn$mean - bysnn$se
bysnn$g <- 'SANDPUMA2 (pplacer)'

pcfull <- read.table("~/git/nrps2/benchmarks/jackknife_fullsetsmile/pc.scored.csv", header=T, sep=",")
fpc <- pcfull[pcfull["method"]=="forced_prediCAT",]
fpc$bin = as.numeric(as.character(cut(fpc$snn, breaks = seq(0,2.5,.5), labels=seq(0,2,.5))))
fpc[fpc$snn >= 2.5, ]$bin <- 2.5
fpc$bin[is.na(fpc$bin)] <- 0
fpc$call <- as.character(fpc$call)
fpc$call[fpc$call == "Y"] <- 1
fpc$call[fpc$call == "N"] <- 0
fpc$call <- as.numeric(fpc$call)
fby <- fpc %>% filter(call_made=="Y") %>% group_by(bin) %>% summarize(N=length(call), mean=mean(call), sd=sd(call), se=sd/sqrt(N))
fby$se[is.na(bysnn$se)] <- 0
fby$ymax <- bysnn$mean + bysnn$se
fby$ymin <- bysnn$mean - bysnn$se
fby$g <- 'SANDPUMA1 (FastTree)'

bysnn.all <- rbind(bysnn, fby)

secolor = "#0C5FE8"

ggplot(bysnn.all, aes(x=bin, y=mean, group=g)) + 
  geom_ribbon(alpha=0.3, aes(ymin=ymin, ymax=ymax, fill=g, color=g), linetype=3, size=0.5) + 
  geom_line(size=1.25, aes(color=g)) + 
  theme_classic() + xlab("SNN Score") + ylab("Fraction Correct") +
  theme(
    axis.title.x = element_text(face="bold", size=14), 
    axis.title.y = element_text(face="bold", size=14), 
    axis.text.x=element_text(color="black", size=12), 
    axis.text.y=element_text(color="black", size=12), 
    plot.title = element_text(face="bold", size=16),
    legend.title = element_blank(),
    legend.position = 'top'
  ) +
  scale_fill_brewer(palette='Dark2') + scale_color_brewer(palette='Dark2')




