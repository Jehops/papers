setwd("/home/jrm/scm/sba.git/doc/paper/data")
rm(list=ls())

## Globals
N <- 50000
##

## Functions

## Build the fp and tp vectors
FpTp <- function(is.ps,fp,tp,nn,np,sp) {
    fpl <- fp
    tpl <- tp

    for (i in 2:length(is.ps)+1) {
        if(!is.ps[i-1]) { fpl[i] <- fpl[i-1]+1; tpl[i] <- tpl[i-1]   }
        if(is.ps[i-1])  { fpl[i] <- fpl[i-1];   tpl[i] <- tpl[i-1]+1 }
    }

    fpl <- fpl/nn
    tpl <- tpl/np

    for (i in 1:sp) {
        fpl[i] <- fpl[sp]
        tpl[i] <- tpl[sp]
    }
    
    eval.parent(substitute(fp<-fpl))
    eval.parent(substitute(tp<-tpl))
}
##

GetPSVec <- function(prefix) {
    ps.sites.file <- paste0("./",prefix,"_ps_sites.csv")
    max.ps.sites <- max(count.fields(ps.sites.file,sep=','))
    ps.sites <- read.table(ps.sites.file,header=F,sep=',',
                           col.names=paste0("V",seq_len(max.ps.sites)),
                           fill=T)
    is.ps <- rep(0,N)
    for (i in 1:100) {
        sim <- ps.sites[i,][!is.na(ps.sites[i,])]
        is.ps[sim + 500*(i-1)] <- 1
    }

    return(is.ps)
}

OrderPS <- function(prefix,model,beb,neb,opt,sba,bebs,nebs,opts,sbas) {
    bebl <- beb; nebl <- neb; optl <- opt; sbal <- sba
    bebsl <- bebs; nebsl <- nebs; optsl <- opts; sbasl <- sbas

    beb.pp <- scan(paste0('./',prefix,'_',model,'_beb_lrt_pp.csv'),sep=',')
    neb.pp <- scan(paste0('./',prefix,'_',model,'_neb_lrt_pp.csv'),sep=',')
    opt.pp <- scan(paste0('./',prefix,'_opt_pp.csv'),sep=',')
    sba.pp <- scan(paste0('./',prefix,'_',model,'_sba_lrt_pp.csv'),sep=',')

    bebl <- bebl[order(1-beb.pp)]; eval.parent(substitute(beb<-bebl))
    nebl <- nebl[order(1-neb.pp)]; eval.parent(substitute(neb<-nebl))
    optl <- optl[order(1-opt.pp)]; eval.parent(substitute(opt<-optl))
    sbal <- sbal[order(1-sba.pp)]; eval.parent(substitute(sba<-sbal))

    bebsl <- length(beb.pp[beb.pp == 1])+1; eval.parent(substitute(bebs<-bebsl))
    nebsl <- length(neb.pp[neb.pp == 1])+1; eval.parent(substitute(nebs<-nebsl))
    optsl <- length(opt.pp[opt.pp == 1])+1; eval.parent(substitute(opts<-optsl))
    sbasl <- length(sba.pp[sba.pp == 1])+1; eval.parent(substitute(sbas<-sbasl))
}

WriteFpTp <- function(prefix,model,fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,
                      tp.opt,tp.sba) {

    n <- length(fp.beb)

    write(fp.beb,paste0('./test/',prefix,'_',model,'_beb_lrt_fp.csv'),ncolumns=n,sep=',')
    write(tp.beb,paste0('./test/',prefix,'_',model,'_beb_lrt_tp.csv'),ncolumns=n,sep=',')
    write(fp.neb,paste0('./test/',prefix,'_',model,'_neb_lrt_fp.csv'),ncolumns=n,sep=',')
    write(tp.neb,paste0('./test/',prefix,'_',model,'_neb_lrt_tp.csv'),ncolumns=n,sep=',')
    write(fp.opt,paste0('./test/',prefix,'_opt_fp.csv'),ncolumns=n,sep=',')
    write(tp.opt,paste0('./test/',prefix,'_opt_tp.csv'),ncolumns=n,sep=',')
    write(fp.sba,paste0('./test/',prefix,'_',model,'_sba_lrt_fp.csv'),ncolumns=n,sep=',')
    write(tp.sba,paste0('./test/',prefix,'_',model,'_sba_lrt_tp.csv'),ncolumns=n,sep=',')
}

### Cat 1, Scheme 3 ###
is.ps <- GetPSVec('c1s3')
N.p <- sum(is.ps)
N.n <- N - N.p

## m2a ##
bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c1s3','m2a',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c1s3','m2a',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

## m8 ##

bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c1s3','m8',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c1s3','m8',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

### Cat 1, Scheme 4 ###
is.ps <- GetPSVec('c1s4')
N.p <- sum(is.ps)
N.n <- N - N.p

## m2a ##
bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c1s4','m2a',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c1s4','m2a',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

## m8 ##

bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c1s4','m8',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c1s4','m8',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

### Cat 2, Scheme 7 ###
is.ps <- GetPSVec('c2s7')
N.p <- sum(is.ps)
N.n <- N - N.p

## m2a ##
bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c2s7','m2a',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c2s7','m2a',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

## m8 ##

bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c2s7','m8',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c2s7','m8',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

### Cat 2, Scheme 8 ###
is.ps <- GetPSVec('c2s8')
N.p <- sum(is.ps)
N.n <- N - N.p

## m2a ##
bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c2s8','m2a',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c2s8','m2a',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

## m8 ##

bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c2s8','m8',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c2s8','m8',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

### Cat 1, Scheme 3, 30 taxa ###
is.ps <- GetPSVec('c1s3_30t')
N.p <- sum(is.ps)
N.n <- N - N.p

## m2a ##
bebo <- nebo <- opto <- sbao <- is.ps
bebs <- nebs <- opts <- sbas <- 1

OrderPS('c1s3_30t','m2a',bebo,nebo,opto,sbao,bebs,nebs,opts,sbas)

fp.beb <- fp.neb <- fp.opt <- fp.sba <- rep(0,N+1)
tp.beb <- tp.neb <- tp.opt <- tp.sba <- rep(0,N+1)

FpTp(bebo,fp.beb,tp.beb,N.n,N.p,bebs)
FpTp(nebo,fp.neb,tp.neb,N.n,N.p,nebs)
FpTp(opto,fp.opt,tp.opt,N.n,N.p,opts)
FpTp(sbao,fp.sba,tp.sba,N.n,N.p,sbas)

WriteFpTp('c1s3_30t','m2a',fp.beb,fp.neb,fp.opt,fp.sba,tp.beb,tp.neb,tp.opt,tp.sba)

#### old plots ####

## Add points to a ROC curve at equally spaced false positive values
## AddPoints <- function(fp,tp,npoints,pch,col) {
##     x <- seq(min(fp),max(fp),length=npoints)
##     id <- rep(npoints,npoints)
##     for(i in 1:(npoints))
##         id[i] <- min(c(1:length(fp))[fp>=x[i]])
##     points(fp[id],tp[id],pch=pch,col=col)
## }

##pdf("~/c1_s3_m2a_new_roc.pdf")

##plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="50% w=1 50% W=1.5 - M2a")
##AddPoints(fp.beb,tp.beb,20,15,'red')
##lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
##AddPoints(fp.neb,tp.neb,20,18,'purple')
##lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
##AddPoints(fp.opt,tp.opt,20,17,'blue')
##lines(fp.sba,tp.sba,type='l',col='green',pch=16)
##AddPoints(fp.sba,tp.sba,20,16,'green')
##abline(0,1)
##legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##pch=c(15,18,17,16))

## dev.off()

##pdf("~/c1_s3_m8_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="50% w=1 50% W=1.5 - M8")
## AddPoints(fp.beb,tp.beb,20,15,'red')
## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')
## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')
## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')
## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
## pch=c(15,18,17,16))

## dev.off()

##pdf("~/c1_s4_m2a_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="45% w=0 45% w=1 10% w=5 - M2a")
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##        pch=c(15,18,17,16))

## dev.off()

## pdf("~/c1_s4_m8_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="45% w=0 45% w=1 10% w=5 - M8")
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##        pch=c(15,18,17,16))

## dev.off()

## pdf("~/c2_s7_m2a_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="50% w=1 50% W=1.5 - M2a - Mild Misspec.")
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##        pch=c(15,18,17,16))

## dev.off()

## pdf("~/c2_s7_m8_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="50% w=1 50% W=1.5 - M8 - Mild Misspec.")
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##        pch=c(15,18,17,16))

## dev.off()

## pdf("~/c2_s8_m2a_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="45% w=0 45% w=1 10% w=5 - M2 - Mild Misspec.",xlim=c(0,.1),ylim=c(0,.8))
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## ## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## ## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## ## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##         pch=c(15,18,17,16))

## dev.off()

## pdf("~/c2_s8_m8_new_roc.pdf")

## plot(fp.beb,tp.beb,type='l',col='red',pch=15,xlab="False Positive Rate",ylab="True Positive Rate",main="45% w=0 45% w=1 10% w=5 - M8 - Mild Misspec.")
## AddPoints(fp.beb,tp.beb,20,15,'red')

## lines(fp.neb,tp.neb,type='l',col='purple',pch=18)
## AddPoints(fp.neb,tp.neb,20,18,'purple')

## lines(fp.opt,tp.opt,type='l',col='blue',pch=17)
## AddPoints(fp.opt,tp.opt,20,17,'blue')

## lines(fp.sba,tp.sba,type='l',col='green',pch=16)
## AddPoints(fp.sba,tp.sba,20,16,'green')

## abline(0,1)
## legend(0.8,0.4,c('BEB','NEB','OPT','SBA'),col=c('red','purple','blue','green'),
##        pch=c(15,18,17,16))

## dev.off()