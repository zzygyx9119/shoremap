#Metod 3... go back to method 1 and add more constraints
library(bbmle)

args<-commandArgs(trailingOnly=TRUE)

high_file<-args[1]
low_file<-args[2]
qtl_file<-args[3]
outprefix<-args[4]
ls.t<-read.table(paste(outprefix,".low.AF.txt",sep=""))
hs.t<-read.table(paste(outprefix,".high.AF.txt",sep=""))
#Parameters
minCoverage <- 4
maxCoverage <- 300
chr<-3

#prep data
hs<-read.table(high_file)
ls<-read.table(low_file)

#extract instersection of markers
chrmod<-10**ceiling(log10(max(hs[,2],ls[,2])))
hs.mod<-hs[,1]*chrmod+hs[,2]
ls.mod<-ls[,1]*chrmod+ls[,2]

intersect.mod<-intersect(hs.mod,ls.mod)
hs<-hs[hs.mod %in% intersect.mod,]
ls<-ls[ls.mod %in% intersect.mod,]


hs.cov<-rowSums(hs[,3:5])
ls.cov<-rowSums(ls[,3:5])

goodCov<-(hs.cov>minCoverage & hs.cov<maxCoverage)&(ls.cov>minCoverage &ls.cov<maxCoverage)

hs<-hs[goodCov,]
ls<-ls[goodCov,]

toUse<-hs[,1]==chr
data<-cbind(ls[toUse,2:4],hs[toUse,3:4])



toUse<-ls.t[,1]==chr

#a.t should contain windowed estimates
a.t<-data.frame(pos=ls.t[toUse,2],ls=ls.t[toUse,3]/100,hs=hs.t[toUse,3]/100)

#use smoothing with zero memmory to get estimates at marker positions
intersect.pos<-intersect(a.t[,1],data[,1])
a.t<-a.t[a.t[,1] %in% intersect.pos,]
data<-data[data[,1] %in% intersect.pos,]

#read qtl file
if(qtl_file==""){
 qtls<-data.frame(pos=-1,chr=-1)
}else{
 qtls<-read.table(qtl_file)
 rank<-sort(sort(abs(qtls[,5]),index.return=TRUE,decreasing=TRUE)$ix,index.return=TRUE)$ix
 lg<-sapply(qtls[,3],function(x) ifelse(sum(qtls[,3]==x)>1,x,0))
 lgType<-sapply(lg,function(x) ifelse(x==0,0,{
  tt<-table(sign(qtls[lg==x,5]))
  ifelse(length(tt)==1,sign(sum(qtls[lg==x,5])),ifelse(tt[1]<tt[2],-tt[1]/tt[2],ifelse(tt[1]==tt[2],0,tt[2]/tt[1])))
 }))
 qtls<-data.frame(id=qtls[,1],pos=qtls[,4],chr=qtls[,3],effect=qtls[,5],rank=rank,lg=lg,lgType=lgType)
}
qtls<-qtls[qtls[,3]==chr,]



singleQtl<-function(q){
# x<-a.t[,1]-q[1]
 q[2]*exp(-4/3*(a.t[,1]-q[1])/ifelse(a.t[,1]<q[1],-q[3],q[4]))
}

freq<-function(x){
 if(length(x)==1){
  #no qtl
  rep(x[1],length(a.t[,1]))
 }else{
  baseline<-x[1]
  qtlest<-matrix(x[-1],byrow=TRUE,ncol=4)
  baseline+rowSums(apply(qtlest,1,singleQtl))
 }
}

separateParameters<-function(x){
 base<-matrix(x[1:2],nrow=1)
 qtlp<-matrix(x[-(1:2)],nrow=7)
 x1<-c(base[1],qtlp[c(1,2:4),])
 x2<-c(base[2],qtlp[c(1,5:7),])
 data.frame(x1=x1,x2=x2)
}

estPlot<-function(x){
 plot(a.t[,c(1,3)],type="l",ylim=c(0,1))
 lines(a.t[,1:2])
 apply(matrix(unique(qtls[,2]),ncol=1),1,function(x) abline(v=x))
 par<-separateParameters(x)
 lines(a.t[,1],freq(par$x1),col="red")
 lines(a.t[,1],freq(par$x2),col="red")
 if(nrow(par)>1){
  apply(matrix(unique(x[seq(3,length(x),7)]),ncol=1),1,function(x) abline(v=x,col="red"))
 }
}

checkRecombinationCount<-function(x,chrsize,verbose=FALSE){
 qtlest<-matrix(x[-1],ncol=4,byrow=TRUE)
 totRec<-sum(qtlest[,1]/qtlest[,3]+(chrsize-qtlest[1])/qtlest[,4])
 if(totRec>100){
  TRUE
 }else{
  FALSE
 }
}

checkRecombinationRates<-function(x,chrsize,verbose=FALSE){
 qtlest<-matrix(x[-1],ncol=4,byrow=TRUE)
 epsilon<-1e-10 #error tolerance...
 if(nrow(qtlest)>1){
  qtlest<-qtlest[sort(qtlest[,1],index.return=TRUE)$ix,]
  #find transitions between negative and positive qtl loci
  transitions<-diff(sign(qtlest[,2]))
  pm<-which(transitions!=0)
  chk<--1
  #demand that the number of recombinations towards the differently directed qtl are at least as many as on the other side
  if(length(pm)>0){
   for(i in pm){
    if(qtlest[i,1]/qtlest[i,3]>(chrsize-qtlest[i,1])/qtlest[i,4]+epsilon){
     if(verbose) print(paste("qtl at position:",qtlest[i,1],"has problems with the recombination rates"))
     chk<-qtlest[i,1]
     break
    }
   }
   if(chk<0){
    for(i in pm+1){
     if(qtlest[i,1]/qtlest[i,3]+epsilon<(chrsize-qtlest[i,1])/qtlest[i,4]){
      if(verbose) print(paste("qtl at position:",qtlest[i,1],"has problems with the recombination rates"))
      chk<-qtlest[i,1]
      break
     }
    }
   }
  }
  chk
 }else{
  #OK
  -1
 }
}

parameterConstraints_sub<-function(df,verbose=FALSE){
 chrsize<-max(a.t[,1])+1
 if(any(abs(df[1,]-0.5)>0.5)){
  if(verbose) print("baseline out of bounds")
  2
 }else if(nrow(df)==1){
  1
 }else if(any(df[seq(2,nrow(df),4),]<0)){
  if(verbose) print("theta too small")
  3
 }else if(any(df[seq(2,nrow(df),4),]>chrsize)){
  if(verbose) print("theta too large")
  4
 }else if(any(abs(df[seq(3,nrow(df),4),])>1)){
  if(verbose) print("effect out of bounds")
  5
 }else if(any(df[c(seq(4,nrow(df),4),seq(5,nrow(df),4)),]<=5e6)){
  if(verbose) print("too small rates")
  6
 }else if(any(abs(df[seq(3,nrow(df),4),1]+df[1,1]-0.5)>0.5)){
  if(verbose) print("first estimate out of bounds")
  7
 }else if(any(abs(df[seq(3,nrow(df),4),2]+df[1,2]-0.5)>0.5)){
  if(verbose) print("second estimate out of bounds")
  8
 }else if(checkRecombinationCount(df$x1,chrsize,verbose)){
  #Failed transition between positive and negative qtls for the first case
  if(verbose) print("first")
  9
 }else if(checkRecombinationCount(df$x2,chrsize,verbose)){
  #Failed transition between positive and negative qtls for the second case
  if(verbose) print("second")
  10
 }else if(checkRecombinationRates(df$x1,chrsize,verbose)>0){
  #Failed transition between positive and negative qtls for the first case
  if(verbose) print("first")
  11
 }else if(checkRecombinationRates(df$x2,chrsize,verbose)>0){
  #Failed transition between positive and negative qtls for the second case
  if(verbose) print("second")
  12
 }else{
  1
 }
}

parameterConstraints<-function(df,verbose=FALSE){
 if(parameterConstraints_sub(df,verbose)==1){
  1
 }else{
  NA
 }
}

getNewQtl<-function(x,direction=c(-1,1)){
 chrsize<-max(a.t[,1])+1
 posToUse<-if(length(x)==2){
  round(chrsize/2)
 }else{
  df<-separateParameters(x)
  #find section with most deviation per marker
  pos<-matrix(sort(df[seq(2,nrow(df),4),1]),ncol=1)
  windows<-rowSums(apply(pos,1,function(p) ifelse(a.t[,1]<p,0,1)))
  distance<-(a.t[,2]-freq(df$x1))**2+(a.t[,3]-freq(df$x2))**2
  round((diff(c(0,pos,chrsize))/2+c(0,pos))[which.max(tapply(distance,windows,mean))])
 }
 #set direction of effects
 effect<-sign(direction)*0.05
 c(posToUse,effect[1],1e7,1e7,effect[2],1e7,1e7)
}

addQtl_sub<-function(x,direction=c(-1,1)){
 chrsize<-max(a.t[,1])+1
 newQtl<-getNewQtl(x,chrsize,direction)
 par<-c(x[1:2],posToUse,effect[1],1e7,1e7,effect[2],1e7,1e7,x[-(1:2)])
 df<-separateParameters(par)
 pc<-parameterConstraints_sub(df)
 if(pc!=1){
  #bad parameters... try to fix
  if(pc==11 ||pc==12){
   chk<-checkRecombinationRates(df$x1,chrsize)
   while(chk>0){
    i<-which(par==chk)
    par[i+3]<-par[i+2]*(chrsize-par[i])/par[i]
    df<-separateParameters(par)
    chk<-checkRecombinationRates(df$x1,chrsize)
   }
   chk<-checkRecombinationRates(df$x2,chrsize)
   while(chk>0){
    i<-which(par==chk)
    par[i+6]=par[i+5]*(chrsize-par[i])/par[i]
    df<-separateParameters(par)
    chk<-checkRecombinationRates(df$x2,chrsize)
   }
   par
  }else{
   #fail 
   print("failed to auto-add a qtl")
   par
  }
 }else{
  par
 }
}

addQtl<-function(previous){
 par.start<-getNewQtl(previous)
 if(is.na(parameterConstraints(separateParameters(c(previous,par.start))))){
  par.start<-getNewQtl(previous,c(1,-1))
 }
 estPlot(c(previous,par.start))
 opt<-optim(freqBoth_addOne, par=par.start, static=previous, method="Nelder-Mead", control=list(parscale=par.start, maxit=5000))
 estPlot(c(previous,opt$par))
 print(opt$val)
 val<-Inf
 while (opt$val<val-1e-3 && opt$convergence==1) {
  val<-opt$val
  opt<-optim(freqBoth, par=opt$par, static=previous, method="Nelder-Mead", control=list(parscale=opt$par, maxit=5000))
  estPlot(c(previous,opt$par))
  print(opt$val)
 }
 opt
}

estLambda<-function(x,col=2){
 sum((a.t[,col]-freq(x))**2)
}

freqBoth<-function(x){
 if(length(x)==2){
   estLambda(x[1],2)+estLambda(x[2],3)
 }else{
  par<-separateParameters(x)
  if(is.na(parameterConstraints(par))){
   NA
  }else{
   estLambda(par$x1,2)+estLambda(par$x2,3)
  }
 }
}

freqBoth_addOne<-function(x,static=c(0.5,0.5)){
 freqBoth(c(static,x))
}

llBoth<-function(x){
 df<-separateParameters(x) #0.001
 if(is.na(parameterConstraints(df))){
  NA
 }else{
#  f1<-freq(df$x1) #0.008
#  f2<-freq(df$x2) #0.007
#  d1<-cbind(data[,2],rowSums(data[,2:3]),f1) #0.007
#  d2<-cbind(data[,4],rowSums(data[,4:5]),f2) #0.007
#  -sum(apply(d1,1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))-sum(apply(d2,1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE))) #0.556
  -sum(apply(cbind(data[,2],rowSums(data[,2:3]),freq(df$x1)),1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))-sum(apply(cbind(data[,4],rowSums(data[,4:5]),freq(df$x2)),1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))
 }
}

calcFocusScale<-function(x){
 if(length(x)>2){
  x*c(1,1,rep(0.5,7),rep(1,length(x)-9))
 }else{
  x
 }
}

calcbic<-function(x){
 2*llBoth(x)+length(x)*log(nrow(data))
}

findBest<-function(par.start,focus=TRUE){
 system.time({
 estPlot(par.start)
 scale<-par.start
 if(focus){
  scale<-calcFocusScale(par.start)
 }
 opt<-optim(freqBoth, par=par.start, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
 estPlot(opt$par)
 print(opt$val)
 val<-Inf
 while (opt$val<val-1e-3 && opt$convergence==1) {
  val<-opt$val
  scale<-opt$par
  if(focus){
   scale<-calcFocusScale(opt$par)
  }
  opt<-optim(freqBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
  estPlot(opt$par)
  print(opt$val)
 }
 })
 print("Refining...")
 par.start<-opt$par

 system.time({
 scale<-par.start
 if(focus){
  scale<-calcFocusScale(par.start)
 }
 opt<-optim(llBoth,par=par.start, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
 estPlot(opt$par)
 print(opt$val)
 val<-Inf
 while (opt$val<val && opt$convergence==1) {
 #while (opt$convergence==1) {
  val<-opt$val
  scale<-opt$par
  if(focus){
   scale<-calcFocusScale(opt$par)
  }
  opt<-optim(llBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
  estPlot(opt$par)
  print(opt$val)
 }
 })
 opt
}



pars<-list()
pars[[1]]<-findBest(c(0.5,0.5))
newQtl<-getNewQtl(pars[[1]]$par)
pars[[2]]<-findBest(c(pars[[1]]$par,newQtl))
#par.start<-c(0.5,0.5,max(a.t[,1])/2,-0.3,1e7,1e7,0.3,1e7,1e7)
#opt<-findBest(par.start)
#pars[[1]]<-opt
for(i in 3:5){
 newQtl<-addQtl(pars[[i-1]]$par)
 pars[[i]]<-findBest(c(pars[[i-1]]$par,newQtl$par))
# t2<-findBest(addQtl(pars[[i-1]]$par,c(1,-1)))
# pars[[i]]<-if(calcbic(t1$par)>calcbic(t2$par)){
#  t2
# }else{
#  t1
# }
}
































#not needed???

ll0w<-function(p1,p2){
 x<-c(p1,p2)
 llBoth(x)
}

ll1w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
 llBoth(x)
}

ll2w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16)
 llBoth(x)
}

ll3w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23)
 llBoth(x)
}

ll4w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21,p22,p23,p24,p25,p26,p27,p28,p29,p30)
 llBoth(x)
}



mle.start<-c(list(),opt$par)
names(mle.start)<-paste("p",1:length(opt$par),sep="")

system.time(m<-mle2(ll1w, start=mle.start))

system.time(m<-mle2(ll2w, start=mle.start))























