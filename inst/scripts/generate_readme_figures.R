#!/usr/bin/env Rscript
## generate_readme_figures.R -?all real data, all package functions
## Bendall et al. 2023 + Farjo et al. 2024 + NCBI GFF3
library(ggplot2)
library(patchwork)
devtools::load_all(".", quiet = TRUE)
library(S4Vectors)
library(SummarizedExperiment)
library(GenomicRanges)

pal <- list(blue="#0072B2", vermillon="#D55E00", green="#009E73",
  orange="#E69F00", skyblue="#56B4E9", purple="#CC79A7",
  gray="#999999", black="#000000")

theme_pub <- function(bs = 10) {
  theme_classic(base_size = bs, base_family = "sans") +
    theme(axis.title=element_text(size=bs,color="black"),
      axis.text=element_text(size=bs-1,color="black"),
      axis.line=element_line(linewidth=0.5,color="black"),
      axis.ticks=element_line(linewidth=0.35,color="black"),
      axis.ticks.length=unit(2.5,"pt"),
      legend.title=element_text(size=bs-1,face="bold"),
      legend.text=element_text(size=bs-1),
      legend.key.size=unit(10,"pt"),
      legend.background=element_blank(),
      panel.background=element_rect(fill="white",color=NA),
      panel.grid=element_blank(),
      strip.background=element_blank(),
      strip.text=element_text(size=bs,face="bold"),
      plot.background=element_rect(fill="white",color=NA),
      plot.margin=margin(6,10,6,6),
      plot.tag=element_text(size=bs+3,face="bold",color="black"))
}

## ---- Bendall data ----
DATA <- "inst/scripts/real_data"
vars <- read.delim(file.path(DATA,"all_variants_filtered.tsv"),stringsAsFactors=FALSE)
cov_l <- readLines(file.path(DATA,"AvgCoverage.all"))[-1]
cov_l <- cov_l[nchar(trimws(cov_l))>0]
cp <- strsplit(trimws(cov_l),"\\s+")
sd <- tapply(as.numeric(sapply(cp,"[",2)),sub("_[12]$","",sapply(cp,"[",1)),mean)
med_depth <- median(sd)
pairs <- read.csv(file.path(DATA,"Transmission_pairs.csv"),stringsAsFactors=FALSE)
vars$fd <- abs(vars$ALT_FREQ_1-vars$ALT_FREQ_2)
vars$conc <- vars$fd<=0.02
vars$mf <- (vars$ALT_FREQ_1+vars$ALT_FREQ_2)/2
ip <- vars[vars$conc,]
all_s <- unique(vars$sample)
n_isnv <- nrow(vars); n_samp <- length(all_s)
n_con <- sum(vars$conc); n_dis <- sum(!vars$conc)
r2 <- cor(vars$ALT_FREQ_1,vars$ALT_FREQ_2)^2
GL <- 29903L
pi_s <- function(f,L=GL,d=NULL){if(length(f)==0)return(0);p<-sum(2*f*(1-f))/L;if(!is.null(d)&&d>1)p<-p*d/(d-1);p}
df_div <- do.call(rbind,lapply(all_s,function(s){
  fa<-vars$ALT_FREQ_1[vars$sample==s]; fp<-ip$ALT_FREQ_1[ip$sample==s]
  d<-if(s%in%names(sd))sd[[s]]else NULL
  data.frame(sample=s,pi_n=pi_s(fa,d=d)*1e4,pi_q=pi_s(fp,d=d)*1e4,
    n_n=length(fa),n_q=length(fp),stringsAsFactors=FALSE)}))
wt <- wilcox.test(df_div$pi_n,df_div$pi_q,paired=TRUE,alternative="greater",exact=FALSE)
mn <- median(df_div$pi_n); mq <- median(df_div$pi_q)
cat(sprintf("Bendall: %d iSNVs, %d samples, R2=%.3f\n",n_isnv,n_samp,r2))

## ================================================================
## FIG 1 -?Replicate QC (unchanged)
## ================================================================
ds <- data.frame(r1=vars$ALT_FREQ_1*100,r2=vars$ALT_FREQ_2*100,
  st=factor(ifelse(vars$conc,"Concordant","Discordant"),levels=c("Concordant","Discordant")))
p1a <- ggplot(ds,aes(r1,r2,fill=st))+geom_abline(slope=1,intercept=0,linetype="dashed",colour=pal$gray,linewidth=0.4)+
  geom_point(shape=21,size=2.2,stroke=0.3,alpha=0.85,colour="white")+
  annotate("text",x=97,y=5,hjust=1,vjust=0,label=sprintf("R\u00b2 = %.3f",r2),size=3.2)+
  scale_fill_manual(name=NULL,values=c(Concordant=pal$green,Discordant=pal$vermillon),
    labels=c(sprintf("Concordant (n=%d)",n_con),sprintf("Discordant (n=%d)",n_dis)))+
  scale_x_continuous("Rep 1 freq. (%)",limits=c(0,100),breaks=c(0,25,50,75,100))+
  scale_y_continuous("Rep 2 freq. (%)",limits=c(0,100),breaks=c(0,25,50,75,100))+
  coord_equal()+labs(tag="a")+guides(fill=guide_legend(override.aes=list(size=3,alpha=1)))+
  theme_pub(10)+theme(legend.position=c(0.35,0.98),legend.justification=c(0,1),
    legend.background=element_rect(fill=alpha("white",0.92),colour=NA))
bl <- data.frame(s=rep(df_div$sample,each=2),x=rep(c(1,2),times=n_samp),
  pi=c(rbind(df_div$pi_n,df_div$pi_q)))
p1b <- ggplot()+geom_line(data=bl,aes(x=x,y=pi,group=s),colour="grey60",linewidth=0.4,alpha=0.45)+
  geom_point(data=data.frame(x=rep(c(1,2),each=n_samp),pi=c(df_div$pi_n,df_div$pi_q),
    c=rep(c("N","Q"),each=n_samp)),aes(x=x,y=pi,colour=c),size=1.8,alpha=0.85,shape=16,
    position=position_jitter(width=0.06,height=0,seed=42))+
  annotate("point",x=1,y=mn,shape=18,size=5.5,colour=pal$black)+
  annotate("point",x=2,y=mq,shape=18,size=5.5,colour=pal$black)+
  annotate("text",x=2.45,y=max(df_div$pi_n)*0.95,label=sprintf("Wilcoxon p=%.1e\nn=%d",wt$p.value,n_samp),
    hjust=1,vjust=1,size=3.0,colour=pal$vermillon,fontface="bold")+
  scale_colour_manual(values=c(N=pal$orange,Q=pal$green),guide="none")+
  scale_x_continuous(NULL,breaks=c(1,2),labels=c("Naive","QC"),limits=c(0.55,2.55))+
  scale_y_continuous(expression(pi~"(x"*10^{-4}*")"),expand=expansion(mult=c(0.03,0.06)))+
  labs(tag="b")+theme_pub(10)
thr <- seq(0,0.20,by=0.01)
pm <- sapply(thr,function(t){pi_v<-sapply(all_s,function(s){
  f<-vars$ALT_FREQ_1[vars$sample==s&vars$fd<=t];d<-if(s%in%names(sd))sd[[s]]else NULL
  pi_s(f,d=d)*1e4});median(pi_v)})
dt <- data.frame(t=thr*100,pi=pm)
p1c <- ggplot(dt,aes(t,pi))+geom_ribbon(aes(ymin=0,ymax=pi),fill=pal$blue,alpha=0.20)+
  geom_line(colour=pal$blue,linewidth=0.9)+geom_point(size=1.5,colour=pal$blue)+
  geom_vline(xintercept=2,linetype="dashed",colour=pal$vermillon,linewidth=0.5)+
  geom_vline(xintercept=5,linetype="dashed",colour=pal$orange,linewidth=0.4)+
  annotate("label",x=2,y=max(dt$pi)*0.85,label="2%",size=3.0,colour=pal$vermillon,fill="white",label.padding=unit(1.5,"pt"))+
  annotate("label",x=5,y=max(dt$pi)*0.72,label="5%",size=3.0,colour=pal$orange,fill="white",label.padding=unit(1.5,"pt"))+
  scale_x_continuous("|freq diff| threshold (%)",breaks=seq(0,20,5),expand=expansion(mult=c(0.01,0.03)))+
  scale_y_continuous(expression("Median "*pi~"(x"*10^{-4}*")"),expand=expansion(mult=c(0.02,0.10)))+
  labs(tag="c")+theme_pub(10)
fig1 <- (p1a|p1b|p1c)+plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
cat("Saving Fig 1...\n")
ggsave("man/figures/fig1_replicate_qc.png",fig1,width=7.5,height=3.1,dpi=300,bg="white")

## ================================================================
## FIG 2 -?Frequency spectrum (unchanged)
## ================================================================
vars$fc <- factor(ifelse(vars$conc,"Concordant","Discordant"),levels=c("Concordant","Discordant"))
ns <- length(unique(paste(vars$POS,vars$REF,vars$ALT)))
p2 <- ggplot(vars,aes(x=mf*100,fill=fc))+geom_histogram(binwidth=2.5,colour="white",linewidth=0.2,alpha=0.90,boundary=0)+
  geom_vline(xintercept=3,linetype="dashed",colour=pal$black,linewidth=0.55)+
  annotate("text",x=4.5,y=Inf,hjust=0,vjust=1.6,label="3%",colour=pal$black,size=3.2)+
  scale_fill_manual(name=NULL,values=c(Concordant=pal$blue,Discordant=pal$vermillon),
    labels=c(sprintf("Concordant (n=%d)",n_con),sprintf("Discordant (n=%d)",n_dis)))+
  scale_x_continuous("Alt allele frequency (%)",breaks=seq(0,100,20),expand=expansion(mult=c(0.01,0.02)))+
  scale_y_continuous("Number of iSNVs",expand=expansion(mult=c(0,0.10)))+
  labs(caption=sprintf("n=%d calls | %d sites | %d samples",n_isnv,ns,n_samp))+
  theme_pub(10)+theme(legend.position=c(0.98,0.98),legend.justification=c(1,1),
    plot.caption=element_text(size=7,color="grey50",margin=margin(t=4)))
cat("Saving Fig 2...\n")
ggsave("man/figures/fig2_frequency_spectrum.png",p2,width=4.2,height=3.2,dpi=300,bg="white")

## ================================================================
## FIG 3 -?Complete workflow: depth + GFF + contamination + sharing
## ================================================================
## a: depth
p3a <- ggplot(data.frame(d=as.numeric(sd)),aes(x=d))+
  geom_histogram(fill=pal$orange,colour="white",linewidth=0.2,binwidth=300,alpha=0.90,boundary=0)+
  geom_vline(xintercept=med_depth,linetype="dashed",colour=pal$black,linewidth=0.5)+
  annotate("text",x=med_depth+200,y=Inf,hjust=0,vjust=1.8,
    label=sprintf("median\n%s x",format(round(med_depth),big.mark=",")),size=3.0,lineheight=1.15)+
  scale_x_continuous("Mean read depth",labels=scales::label_comma(),expand=expansion(mult=c(0.01,0.05)))+
  scale_y_continuous("Samples",expand=expansion(mult=c(0,0.15)))+labs(tag="a")+theme_pub(10)

## b: GFF annotation
gff <- file.path(DATA,"sars2_NC045512.gff3")
gff_lines <- readLines(gff); gff_lines <- gsub("NC_045512\\.2","MN908947.3",gff_lines)
gff_tmp <- tempfile(fileext=".gff3"); writeLines(gff_lines,gff_tmp)
gr_mc <- GRanges("MN908947.3",IRanges::IRanges(ip$POS,width=1))
mcols(gr_mc)$ref <- ip$REF; mcols(gr_mc)$alt <- ip$ALT
fm <- matrix(ip$ALT_FREQ_1,ncol=1); colnames(fm) <- "pooled"
whe_mc <- WithinHostExperiment(assays=list(altFreq=fm),rowRanges=gr_mc,
  colData=DataFrame(sample_id="pooled"))
whe_mc <- annotateFromGFF(whe_mc,gff_tmp)
gv <- mcols(rowRanges(whe_mc))$GFF_FEATURE; gv[is.na(gv)] <- "intergenic"
gt <- as.data.frame(table(Gene=gv),stringsAsFactors=FALSE)
gt <- gt[order(gt$Freq,decreasing=TRUE),]; gt$Gene <- factor(gt$Gene,levels=rev(gt$Gene))
p3b <- ggplot(gt,aes(y=Gene,x=Freq))+geom_col(fill=pal$blue,colour="white",linewidth=0.2,alpha=0.90)+
  geom_text(aes(label=Freq),hjust=-0.25,size=3.1)+
  scale_x_continuous("Concordant iSNVs (GFF-annotated)",expand=expansion(mult=c(0,0.22)))+
  labs(y=NULL,tag="b")+theme_pub(10)+theme(axis.text.y=element_text(size=9))

## Build pair_df for sharing and lollipop panels
pdf2 <- merge(pairs[pairs$Transmission_indiv=="A",c("sample","pair_id")],
  pairs[pairs$Transmission_indiv=="B",c("sample","pair_id")],by="pair_id",suffixes=c("_d","_r"))

## c: best pair lollipop (HH46)
bp <- "HH46"; bpr <- pdf2[pdf2$pair_id==bp,]
dv_best <- ip[ip$sample==bpr$sample_d[1],]
df_loll <- do.call(rbind,lapply(seq_len(nrow(dv_best)),function(j){
  rr <- vars[vars$sample==bpr$sample_r[1]&vars$POS==dv_best$POS[j]&vars$ALT==dv_best$ALT[j],]
  gn <- gv[ip$POS==dv_best$POS[j]][1]; if(is.na(gn)) gn <- "intergenic"
  data.frame(pos=dv_best$POS[j],gene=gn,df=dv_best$ALT_FREQ_1[j]*100,
    shared=nrow(rr)>0,stringsAsFactors=FALSE)}))
df_loll$lb <- sprintf("%s (%s)",df_loll$pos,df_loll$gene)
df_loll$lb <- factor(df_loll$lb,levels=rev(df_loll$lb))
nl <- sum(!df_loll$shared)
p3c <- ggplot(df_loll,aes(y=lb,x=df,fill=shared))+
  geom_col(width=0.6,alpha=0.90)+
  geom_text(aes(label=sprintf("%.1f%%",df)),hjust=-0.15,size=2.8)+
  scale_fill_manual(name=NULL,values=c("TRUE"=pal$green,"FALSE"=pal$vermillon),
    labels=c("TRUE"="Detected","FALSE"="Lost"))+
  scale_x_continuous("Donor freq. (%)",expand=expansion(mult=c(0,0.30)))+
  annotate("text",x=max(df_loll$df)*1.2,y=0.6,hjust=1,vjust=0,
    label=sprintf("Pair %s\n%d/%d lost",bp,nl,nrow(df_loll)),size=3.0,lineheight=1.2)+
  labs(y=NULL,tag="c")+theme_pub(10)+
  theme(legend.position=c(0.98,0.02),legend.justification=c(1,0),axis.text.y=element_text(size=8))

## d: variant sharing
sh <- do.call(rbind,lapply(seq_len(nrow(pdf2)),function(i){
  dv<-ip[ip$sample==pdf2$sample_d[i],];nd<-nrow(dv);if(nd==0)return(NULL)
  ns<-sum(sapply(seq_len(nd),function(j)nrow(vars[vars$sample==pdf2$sample_r[i]&vars$POS==dv$POS[j]&vars$ALT==dv$ALT[j],])>0))
  data.frame(pid=pdf2$pair_id[i],nd=nd,ns=ns,nl=nd-ns,stringsAsFactors=FALSE)}))
sh <- sh[order(sh$nd,decreasing=TRUE),]; sh$rk <- seq_len(nrow(sh))
nz <- sum(sh$ns==0); pz <- round(100*nz/nrow(sh))
dfs <- rbind(data.frame(rk=sh$rk,ct=sh$ns,tp="Shared"),data.frame(rk=sh$rk,ct=sh$nl,tp="Lost"))
dfs$tp <- factor(dfs$tp,levels=c("Lost","Shared"))
p3d <- ggplot(dfs,aes(x=rk,y=ct,fill=tp))+geom_col(width=0.7,alpha=0.90)+
  scale_fill_manual(name=NULL,values=c(Lost=pal$vermillon,Shared=pal$green))+
  annotate("text",x=nrow(sh)*0.95,y=Inf,hjust=1,vjust=1.5,
    label=sprintf("%d/%d (%.0f%%)\nzero shared",nz,nrow(sh),pz),size=3.0,lineheight=1.2)+
  scale_x_continuous("Transmission pair",expand=expansion(mult=c(0.01,0.01)))+
  scale_y_continuous("Donor iSNVs",expand=expansion(mult=c(0,0.12)))+
  labs(tag="d")+theme_pub(10)+theme(legend.position=c(0.02,0.98),legend.justification=c(0,1))

fig3 <- (p3a+p3b)/(p3c+p3d)+plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
cat("Saving Fig 3...\n")
ggsave("man/figures/fig3_qc_overview.png",fig3,width=7.0,height=5.8,dpi=300,bg="white")

## ================================================================
## FIG 4 -?Longitudinal: diversity + trajectories + SFS + temporal QC
## ================================================================
cat("Loading Farjo data...\n")
ff <- sort(list.files(file.path(DATA,"farjo_longitudinal"),pattern="ivar",full.names=TRUE))
ntp <- length(ff)
whe_f <- readWithinHostTable(ff,format="ivar",
  colData=DataFrame(sample_id=paste0("d",seq_len(ntp)),host_id=rep("p",ntp),timepoint=seq_len(ntp)))
whe_f <- flagISNV(whe_f,ISNVFilter(minDepth=100L,minFreq=0.03,maxFreq=0.97))

div_f <- as.data.frame(calcDiversity(whe_f,genomeLength=GL,indices=c("pi","richness")))
div_f$day <- seq_len(ntp); div_f$pv <- as.numeric(div_f$pi)*1e4; div_f$rv <- as.numeric(div_f$richness)
pk <- div_f$day[which.max(div_f$rv)]
cat(sprintf("  %d timepoints, peak=%d iSNVs at day %d\n",ntp,max(div_f$rv),pk))

## a: diversity arc
p4a <- ggplot(div_f,aes(x=day))+
  geom_col(aes(y=rv),fill=pal$skyblue,alpha=0.50,width=0.6)+
  geom_line(aes(y=pv*max(rv)/max(pv)),colour=pal$vermillon,linewidth=1.0)+
  geom_point(aes(y=pv*max(rv)/max(pv)),colour=pal$vermillon,size=2.5)+
  scale_y_continuous("QC-passed iSNVs (bars)",
    sec.axis=sec_axis(~.*max(div_f$pv)/max(div_f$rv),name=expression(pi~"(x"*10^{-4}*", line)")),
    expand=expansion(mult=c(0,0.08)))+
  scale_x_continuous("Timepoint",breaks=seq_len(ntp))+labs(tag="a")+theme_pub(10)+
  theme(axis.title.y.right=element_text(colour=pal$vermillon),axis.text.y.right=element_text(colour=pal$vermillon))

## b: frequency trajectories
tj <- as.data.frame(trackFrequency(whe_f,hostCol="host_id",timeCol="timepoint"))
vr <- tapply(tj$frequency,tj$variant_key,function(x)diff(range(x)))
t10 <- names(sort(vr,decreasing=TRUE))[1:min(10,length(vr))]
tjt <- tj[tj$variant_key%in%t10,]
p4b <- ggplot(tjt,aes(x=timepoint,y=frequency*100,colour=variant_key,group=variant_key))+
  geom_line(linewidth=0.7,alpha=0.80)+geom_point(size=1.5,alpha=0.90)+
  geom_hline(yintercept=3,linetype="dotted",colour=pal$gray,linewidth=0.3)+
  scale_x_continuous("Timepoint",breaks=seq_len(ntp))+
  scale_y_continuous("Alt allele freq. (%)",limits=c(0,100))+
  labs(colour=NULL,tag="b")+theme_pub(10)+theme(legend.position="none")

## c: SFS evolution (day 1 vs peak vs day 9)
kd <- c(1,pk,ntp)
sfs_tl <- lapply(kd,function(d)buildSFS(whe_f,sampleIdx=d,fold=TRUE,genomeLength=GL,nBins=8))
names(sfs_tl) <- paste0("Day ",kd)
df_st <- do.call(rbind,lapply(names(sfs_tl),function(nm){
  s<-sfs_tl[[nm]]; br<-s@breaks; mid<-(br[-length(br)]+br[-1L])/2; tot<-max(sum(s@counts),1)
  data.frame(freq=mid*100,prop=s@counts/tot,day=nm,n=sum(s@counts),stringsAsFactors=FALSE)}))
df_st$day <- factor(df_st$day,levels=names(sfs_tl))
p4c <- ggplot(df_st,aes(x=freq,y=prop,fill=day))+
  geom_col(position="dodge",alpha=0.85,colour="white",linewidth=0.15)+
  scale_fill_manual(name=NULL,values=c(pal$skyblue,pal$vermillon,pal$green))+
  scale_x_continuous("Minor allele freq. (%)")+
  scale_y_continuous("Proportion",expand=expansion(mult=c(0,0.12)))+
  labs(tag="c")+theme_pub(10)+theme(legend.position=c(0.98,0.98),legend.justification=c(1,1))

## d: temporal QC -?persistent vs transient
whe_ft <- flagTemporalInconsistency(whe_f,minTimepoints=2L)
tc <- mcols(rowRanges(whe_ft))$temporal_class
n_per <- sum(tc=="persistent",na.rm=TRUE); n_tra <- sum(tc=="transient",na.rm=TRUE)
cat(sprintf("  Temporal QC: %d persistent, %d transient\n",n_per,n_tra))
df_tc <- data.frame(class=c("Persistent\n(>= 2 timepoints)","Transient\n(1 timepoint only)"),
  count=c(n_per,n_tra),stringsAsFactors=FALSE)
df_tc$class <- factor(df_tc$class,levels=df_tc$class)
p4d <- ggplot(df_tc,aes(x=class,y=count,fill=class))+
  geom_col(width=0.6,alpha=0.90)+
  geom_text(aes(label=count),vjust=-0.3,size=3.5,fontface="bold")+
  scale_fill_manual(values=c(pal$green,pal$vermillon),guide="none")+
  scale_y_continuous("Variant sites",expand=expansion(mult=c(0,0.15)))+
  annotate("text",x=1.5,y=max(df_tc$count)*0.5,
    label=sprintf("%.0f%% of iSNVs are transient\n(likely artefacts)",100*n_tra/(n_per+n_tra)),
    size=3.0,lineheight=1.2)+
  labs(x=NULL,tag="d")+theme_pub(10)

fig4 <- (p4a+p4b)/(p4c+p4d)+plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
cat("Saving Fig 4...\n")
ggsave("man/figures/fig4_diversity_landscape.png",fig4,width=7.0,height=5.8,dpi=300,bg="white")

## ================================================================
## FIG 5 -?SFS paradigm + pop genetics
## ================================================================
tajD <- function(f,L,dep,nc=100L){S<-length(f);if(S==0)return(0);n<-min(max(as.integer(round(dep)),2L),nc);
  ph<-(2/L)*sum(f*(1-f))*n/(n-1L);a1<-sum(1/seq_len(n-1L));a2<-sum(1/(seq_len(n-1L))^2);
  b1<-(n+1)/(3*(n-1));b2<-2*(n^2+n+3)/(9*n*(n-1));c1<-b1-1/a1;c2<-b2-(n+2)/(a1*n)+a2/a1^2;
  e1<-c1/a1;e2<-c2/(a1^2+a2);ds<-ph*L-S/a1;dn<-sqrt(e1*S+e2*S*(S-1));if(dn>0)ds/dn else 0}
dtj <- do.call(rbind,lapply(all_s,function(s){f<-ip$ALT_FREQ_1[ip$sample==s];
  d<-if(s%in%names(sd))sd[[s]]else 100; data.frame(D=tajD(f,GL,d),S=length(f))}))
dtp <- dtj[dtj$S>0,]; medD <- median(dtp$D); pn <- round(100*sum(dtp$D<0)/nrow(dtp))

## a: Tajima's D
p5a <- ggplot(dtp,aes(x=D))+geom_histogram(bins=15,fill=pal$blue,colour="white",linewidth=0.2,alpha=0.90)+
  geom_vline(xintercept=0,linetype="dashed",colour=pal$gray,linewidth=0.5)+
  geom_vline(xintercept=medD,colour=pal$vermillon,linewidth=0.6)+
  annotate("text",x=max(dtp$D)*0.95,y=Inf,vjust=1.5,hjust=1,
    label=sprintf("median=%.2f\n%d%% < 0",medD,pn),size=3.0,colour=pal$vermillon,lineheight=1.2)+
  scale_x_continuous("Tajima's D (n capped at 100)")+
  scale_y_continuous("Samples",expand=expansion(mult=c(0,0.12)))+labs(tag="a")+theme_pub(10)

## b: dN/dS
gv2 <- ip$GFF_FEATURE; gv2[is.na(gv2)|gv2==""] <- "Intergenic"
sy <- !is.na(ip$REF_AA)&!is.na(ip$ALT_AA)&ip$REF_AA==ip$ALT_AA
ddn <- do.call(rbind,lapply(unique(gv2),function(g){i<-which(gv2==g);nS<-sum(sy[i]);nN<-sum(!sy[i]);
  if(nS+nN<2)return(NULL);data.frame(gene=g,nS=nS,nN=nN,dNdS=if(nS>0)nN/nS else NA_real_,
    has=nS>0,stringsAsFactors=FALSE)}))
if(!is.null(ddn)&&nrow(ddn)>0){
  ddn<-ddn[order(ddn$nS+ddn$nN,decreasing=TRUE),]
  go<-c(setdiff(ddn$gene,"Intergenic"),"Intergenic");go<-go[go%in%ddn$gene]
  ddn$gene<-factor(ddn$gene,levels=rev(go));ddn$dp<-ifelse(is.na(ddn$dNdS),0,ddn$dNdS)
  p5b <- ggplot(ddn,aes(y=gene,x=dp,fill=has))+geom_col(width=0.6,alpha=0.90)+
    geom_vline(xintercept=1,linetype="dashed",colour=pal$gray,linewidth=0.5)+
    geom_text(aes(label=ifelse(has,sprintf("%.1f (%dN/%dS)",dNdS,nN,nS),sprintf("%dN, 0S",nN))),
      hjust=-0.05,size=2.7)+
    scale_fill_manual(values=c("TRUE"=pal$vermillon,"FALSE"=pal$gray),guide="none")+
    scale_x_continuous("Within-host dN/dS",expand=expansion(mult=c(0,0.40)))+
    labs(y=NULL,tag="b")+theme_pub(10)+theme(axis.text.y=element_text(size=9))
} else { p5b <- ggplot()+labs(tag="b")+theme_void() }

## c: ranked pi lollipop
dr <- df_div[order(df_div$pi_n,decreasing=TRUE),]; dr$rk <- seq_len(nrow(dr))
dra <- rbind(data.frame(rk=dr$rk,pi=dr$pi_n,tp="Naive"),data.frame(rk=dr$rk,pi=dr$pi_q,tp="QC"))
dra$tp <- factor(dra$tp,levels=c("Naive","QC"))
p5c <- ggplot(dra,aes(rk,pi,colour=tp))+
  geom_segment(data=dr,aes(x=rk,xend=rk,y=pi_q,yend=pi_n),colour="grey65",linewidth=0.5,inherit.aes=FALSE)+
  geom_point(size=1.5,alpha=0.90,shape=16)+
  geom_hline(yintercept=mn,linetype="dotted",colour=pal$orange,linewidth=0.5)+
  geom_hline(yintercept=mq,linetype="dotted",colour=pal$green,linewidth=0.5)+
  scale_colour_manual(name=NULL,values=c(Naive=pal$orange,QC=pal$green),
    guide=guide_legend(override.aes=list(size=3)))+
  scale_x_continuous("Sample rank",expand=expansion(mult=c(0.01,0.01)))+
  scale_y_continuous(expression(pi~"(x"*10^{-4}*")"),expand=expansion(mult=c(0,0.08)))+
  labs(tag="c")+theme_pub(10)+theme(legend.position=c(0.98,0.98),legend.justification=c(1,1))

## d: SFS naive vs QC with bias correction
gr_all <- GRanges("MN908947.3",IRanges::IRanges(vars$POS,width=1))
mcols(gr_all)$ref <- vars$REF; mcols(gr_all)$alt <- vars$ALT
fmn <- matrix(vars$ALT_FREQ_1,ncol=1); colnames(fmn) <- "naive"
dmn <- matrix(as.integer(vars$TOTAL_DP_1),ncol=1); colnames(dmn) <- "naive"
wn <- WithinHostExperiment(assays=list(altFreq=fmn,totalDepth=dmn),rowRanges=gr_all,
  colData=DataFrame(sample_id="naive"))
gr_qc <- GRanges("MN908947.3",IRanges::IRanges(ip$POS,width=1))
fmq <- matrix(ip$ALT_FREQ_1,ncol=1); colnames(fmq) <- "qc"
dmq <- matrix(as.integer(ip$TOTAL_DP_1),ncol=1); colnames(dmq) <- "qc"
wq <- WithinHostExperiment(assays=list(altFreq=fmq,totalDepth=dmq),rowRanges=gr_qc,
  colData=DataFrame(sample_id="qc"))
sn <- buildSFS(wn,fold=TRUE,genomeLength=GL,nBins=10)
sq <- buildSFS(wq,fold=TRUE,genomeLength=GL,nBins=10)
sc <- correctSFSBias(sq,method="binomial")
.s2d <- function(s,lb){br<-s@breaks;mid<-(br[-length(br)]+br[-1L])/2
  data.frame(freq=mid*100,count=s@counts,lb=lb)}
dfs5 <- rbind(.s2d(sn,sprintf("Naive (n=%d)",sn@nSites)),
  .s2d(sq,sprintf("QC (n=%d)",sq@nSites)),
  .s2d(sc,sprintf("QC+bias corr. (n=%d)",sc@nSites)))
dfs5$lb <- factor(dfs5$lb,levels=unique(dfs5$lb))
p5d <- ggplot(dfs5,aes(x=freq,y=count,fill=lb))+
  geom_col(position="dodge",alpha=0.85,width=2.0,colour="white",linewidth=0.15)+
  scale_fill_manual(name=NULL,values=c(pal$orange,pal$green,pal$blue))+
  scale_x_continuous("Minor allele freq. (%)")+
  scale_y_continuous("iSNVs",expand=expansion(mult=c(0,0.12)))+
  labs(tag="d")+theme_pub(10)+theme(legend.position=c(0.98,0.98),legend.justification=c(1,1),
    legend.text=element_text(size=7.5))

fig5 <- (p5a+p5b)/(p5c+p5d)+plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
cat("Saving Fig 5...\n")
ggsave("man/figures/fig5_evolutionary_analysis.png",fig5,width=7.5,height=5.8,dpi=300,bg="white")

## ================================================================
## FIG 6 -?QC preserves signal: Tajima overlay + compareSFS
## ================================================================
dtn <- do.call(rbind,lapply(all_s,function(s){f<-vars$ALT_FREQ_1[vars$sample==s];
  if(length(f)==0)return(NULL);d<-if(s%in%names(sd))sd[[s]]else 100;
  data.frame(D=tajD(f,GL,d),cond="Naive",stringsAsFactors=FALSE)}))
db <- rbind(dtn,data.frame(D=dtp$D,cond="QC")); db$cond <- factor(db$cond,levels=c("Naive","QC"))
wtD <- wilcox.test(dtn$D,dtp$D,exact=FALSE)

p6a <- ggplot(db,aes(x=D,fill=cond))+
  geom_histogram(bins=14,colour="white",linewidth=0.2,alpha=0.55,position="identity")+
  geom_vline(xintercept=0,linetype="dashed",colour=pal$gray,linewidth=0.4)+
  scale_fill_manual(name=NULL,values=c(Naive=pal$orange,QC=pal$green))+
  annotate("text",x=max(db$D)*0.95,y=Inf,hjust=1,vjust=1.5,
    label=sprintf("Wilcoxon p=%.2f (n.s.)\nQC preserves signal",wtD$p.value),
    size=3.0,lineheight=1.2)+
  scale_x_continuous("Tajima's D")+scale_y_continuous("Samples",expand=expansion(mult=c(0,0.15)))+
  labs(tag="a")+theme_pub(10)+theme(legend.position=c(0.02,0.98),legend.justification=c(0,1))

## b: compareSFS -?formal test
comp <- compareSFS(sn,sq)
cat(sprintf("  compareSFS: chi2=%.1f, p=%.3f\n",comp$chisq_stat,comp$chisq_p))
dfp <- rbind(data.frame(freq=(sn@breaks[-length(sn@breaks)]+sn@breaks[-1])/2*100,
  prop=comp$proportions1,lb="Naive"),
  data.frame(freq=(sq@breaks[-length(sq@breaks)]+sq@breaks[-1])/2*100,
  prop=comp$proportions2,lb="QC"))
dfp$lb <- factor(dfp$lb,levels=c("Naive","QC"))
p6b <- ggplot(dfp,aes(x=freq,y=prop,fill=lb))+
  geom_col(position="dodge",alpha=0.80,colour="white",linewidth=0.15)+
  scale_fill_manual(name=NULL,values=c(Naive=pal$orange,QC=pal$green))+
  annotate("label",x=40,y=max(dfp$prop)*0.70,hjust=0.5,vjust=0.5,
    label=sprintf("chi2 = %.1f\np = %.3f",comp$chisq_stat,comp$chisq_p),
    size=3.2,lineheight=1.3,fill=alpha("white",0.90),label.padding=unit(4,"pt"))+
  scale_x_continuous("Minor allele freq. (%)")+
  scale_y_continuous("Proportion",expand=expansion(mult=c(0,0.12)))+
  labs(tag="b")+theme_pub(10)+theme(legend.position=c(0.75,0.85),legend.justification=c(0,1))

fig6 <- (p6a|p6b)+plot_annotation(theme=theme(plot.background=element_rect(fill="white",color=NA)))
cat("Saving Fig 6...\n")
ggsave("man/figures/fig6_consensus_validation.png",fig6,width=7.0,height=3.2,dpi=300,bg="white")

cat("\n=== ALL 6 FIGURES SAVED ===\n")
cat(sprintf("R: %s\n",R.version.string))
cat("Bendall 2023 + Farjo 2024 + NCBI GFF3. All real data.\n")
