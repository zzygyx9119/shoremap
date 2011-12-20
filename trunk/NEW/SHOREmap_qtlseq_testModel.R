#Metod 3... go back to method 1 and add more constraints
#wrapper v01
source("/projects/dep_coupland/grp_nordstrom/projects/Lotus/Simulation/PopSimulatorRIL/Rcode/SHOREmap_qtlseq_testModel_lib.R")

args<-commandArgs(trailingOnly=TRUE)

high_file<-args[1]
low_file<-args[2]
qtl_file<-args[3]
outprefix<-args[4]
#ls.t<-read.table(paste(outprefix,".low.AF.txt",sep=""))
#hs.t<-read.table(paste(outprefix,".high.AF.txt",sep=""))
#Parameters
minCoverage <- 4
maxCoverage <- 300
winSize <- 2e6
tolerance <- 2e6


#prep data
hs_all<-read.table(high_file)
ls_all<-read.table(low_file)

#extract instersection of markers
chrmod<-10**ceiling(log10(max(hs_all[,2],ls_all[,2])))
hs.mod<-hs_all[,1]*chrmod+hs_all[,2]
ls.mod<-ls_all[,1]*chrmod+ls_all[,2]

intersect.mod<-intersect(hs.mod,ls.mod)
hs_all<-hs_all[hs.mod %in% intersect.mod,]
ls_all<-ls_all[ls.mod %in% intersect.mod,]


hs.cov<-rowSums(hs_all[,3:5])
ls.cov<-rowSums(ls_all[,3:5])

goodCov<-(hs.cov>minCoverage & hs.cov<maxCoverage)&(ls.cov>minCoverage &ls.cov<maxCoverage)

hs_all<-hs_all[goodCov,]
ls_all<-ls_all[goodCov,]

#read qtl file
qtls_all<-data.frame(pos=-1,chr=-1)
if(qtl_file==""){
 qtls_all<-data.frame(pos=-1,chr=-1)
}else{
 qtls_all<-read.table(paste(qtl_file,".qtl.txt",sep=""))
 ls.t<-read.table(paste(qtl_file,".low.AF.txt",sep=""))
 hs.t<-read.table(paste(qtl_file,".high.AF.txt",sep=""))
 a.t<-cbind(ls.t[,1],ls.t[,2],(hs.t[,3]-ls.t[,3])/100)
 marker<-qtls_all[,2]+1
 freqDiff<-a.t[marker,3]
 trueStart<-sapply(1:length(freqDiff),function(i){
  start<-marker[i]
  while(a.t[start,3]==freqDiff[i]&&start>1){
   start<-start-1
  }
  a.t[start,2]+1
 })
 trueEnd<-sapply(1:length(freqDiff),function(i){
  end<-marker[i]
  while(a.t[end,3]==freqDiff[i]&& end<nrow(a.t)){
   end<-end+1
  }
  a.t[end,2]-1
 })
 rank<-sort(sort(abs(qtls_all[,5]),index.return=TRUE,decreasing=TRUE)$ix,index.return=TRUE)$ix
 lg<-sapply(qtls_all[,3],function(x) ifelse(sum(qtls_all[,3]==x)>1,x,0))
 lgType<-sapply(lg,function(x) ifelse(x==0,0,{
  tt<-table(sign(qtls_all[lg==x,5]))
  ifelse(length(tt)==1,sign(sum(qtls_all[lg==x,5])),ifelse(tt[1]<tt[2],-tt[1]/tt[2],ifelse(tt[1]==tt[2],0,tt[2]/tt[1])))
 }))
 qtls_all<-data.frame(id=qtls_all[,1],pos=qtls_all[,4],chr=qtls_all[,3],effect=qtls_all[,5],rank=rank,lg=lg,lgType=lgType,trueFreqDiff=freqDiff,trueStart=trueStart,trueEnd=trueEnd)
 
}

verbose<-FALSE

estimates<- matrix(ncol=5)[FALSE,]
colnames(estimates)<-c("chr","start","stop","freqDiff","Est")

nrOfChrs<-length(unique(hs_all[,1]))

pdf(paste(outprefix,".plots.pdf",sep=""))
par(mfrow=c(ceiling(nrOfChrs/2),2))


for (chr in unique(hs_all[,1])){
 print(paste("Chr:",chr))
 if(verbose) {
  X11()
 }
 toUse<-hs_all[,1]==chr
 data<-cbind(ls_all[toUse,2:4],hs_all[toUse,3:4])
 #smoothing... not needed???
 a.t<-t(sapply(data[,1],function(x) {index<-which(data[,1]>=x-winSize/2 & data[,1]<x+winSize/2);c(x,sum(data[index,2])/sum(data[index,2:3]),sum(data[index,4])/sum(data[index,4:5]))}))
 #begin with frequencies... later implement probit??
 #a.t.raw<-cbind(x=data[,1],y1=data[,2]/rowSums(data[,2:3]),y2=data[,4]/rowSums(data[,4:5]))
# using true values
# toUse<-ls.t[,1]==chr
# a.t<-data.frame(pos=ls.t[toUse,2],ls=ls.t[toUse,3]/100,hs=hs.t[toUse,3]/100)
 #use smoothing with zero memmory to get estimates at marker positions
# intersect.pos<-intersect(a.t[,1],data[,1])
# a.t<-a.t[a.t[,1] %in% intersect.pos,]
# data<-data[data[,1] %in% intersect.pos,]


 qtls<-qtls_all[qtls_all[,3]==chr,]
  

 pars<-list()
 pars[[1]]<-findBest(c(0.5,0.5),verbose=verbose)
 if(verbose) print(paste("No QTL:",calcbic(pars[[1]]$par)))
 plus<-findBest(c(pars[[1]]$par,getNewQtl(pars[[1]]$par,c(-1,1))),verbose=verbose)
 minus<-findBest(c(pars[[1]]$par,getNewQtl(pars[[1]]$par,c(1,-1))),verbose=verbose)
 if(calcbic(plus$par)<calcbic(minus$par)){
  pars[[2]]<-plus
 }else{
  pars[[2]]<-minus
 }
 if(verbose) print(paste("1 QTL:",calcbic(pars[[2]]$par)))
 if(calcbic(pars[[2]]$par)<calcbic(pars[[1]]$par)){
  #add QTLs as long as BIC improves
  bestVal<-calcbic(pars[[2]]$par)
  for(i in 3:10){
   plus<-findBest(c(pars[[i-1]]$par,getNewQtl(pars[[i-1]]$par,c(-1,1))),verbose=verbose)
   minus<-findBest(c(pars[[i-1]]$par,getNewQtl(pars[[i-1]]$par,c(1,-1))),verbose=verbose)
   if(calcbic(plus$par)<calcbic(minus$par)){
    pars[[i]]<-plus
   }else{
    pars[[i]]<-minus
   } 
#   newQtl<-getNewQtl(pars[[i-1]]$par,c(-1,1))
#   pars[[i]]<-findBest(c(pars[[i-1]]$par,newQtl),verbose=verbose)
   if(verbose) print(paste(i-1,"QTLs:",calcbic(pars[[i]]$par)))
   newVal<-calcbic(pars[[i]]$par)
   if(newVal<=bestVal){
    bestVal<-newVal
   }else{
    break
   }
  }
 }
 bestPar<-pars[[which.min(sapply(pars,function(x) calcbic(x$par)))]]$par 
 #do plot
 if(verbose) dev.off()
 estPlot(bestPar)
 #output parameters to file
 sapply(pars,function(x) write.table( matrix(c(chr,calcbic(x$par),x$par),nrow=1),file=paste(outprefix,"_modelParameters.csv",sep=""), append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE))
 #add estimates estimate list
 if(length(bestPar)>2){
  estimates<-rbind(estimates,t(sapply(bestPar[seq(3,length(bestPar),7)],function(x) {
   closestMarker<-which.min(abs(a.t[,1]-x))
   c(chr,max(x-tolerance,0),min(x+tolerance,max(a.t[,1])),abs(a.t[closestMarker,2]-a.t[closestMarker,3]),x)
  })))
 }
}

dev.off()


#print estimates
qtls<-qtls_all
header<-c("chr",colnames(qtls)[c(1:2,4:ncol(qtls))],colnames(estimates)[c(2:ncol(estimates))],"judgement","spec")
header[1]<-paste("#",header[1],sep="")
qhc<-ncol(qtls)-1
ehc<-ncol(estimates)-1

toPrint<-sapply(unique(hs_all[,1]),function(chr){
 cq<-sum(qtls$chr==chr)
 ce<-sum(estimates[,1]==chr)
 q<-matrix(c(qtls[qtls$chr==chr,c(1:2,4:ncol(qtls))],recursive=TRUE),ncol=qhc)
 e<-matrix(estimates[estimates[,1]==chr,c(2:ncol(estimates))],ncol=ehc)
 if(cq>0 && ce>0){
  #get True Positive
  pairs<-combn(cq+ce,2)
  pairs<-matrix(pairs[,pairs[1,]<=cq & pairs[2,]>cq],nrow=2)
  pairs[2,]<-pairs[2,]-cq
  tp<-apply(pairs,2,function(x){
   q[x[1],2]>=e[x[2],1] && q[x[1],2]<=e[x[2],2]
  })
  toPrint<-matrix(apply(matrix(pairs[,tp],nrow=2),2,function(x) c(chr,q[x[1],],e[x[2],],"TP","")),byrow=TRUE,ncol=length(header))
  #check if any of the qtls are close to a interval
  close<-apply(pairs,2,function(x){
   q[x[1],2]+tolerance>=e[x[2],1] && q[x[1],2]-tolerance<=e[x[2],2]
  })
  toPrint<-rbind(toPrint,matrix(apply(matrix(pairs[,close &!tp],nrow=2),2,function(x) c(chr,q[x[1],],e[x[2],],"FP","close")),byrow=T,ncol=length(header)))
  #more FP?
  for(i in 1:ce){
   if(!i %in% unique(pairs[2,close]) ){
    toPrint<-rbind(toPrint,matrix(c(chr,rep(NA,qhc),e[i,],"FP","" ),ncol=length(header)))
   }
  }
  #FN
  for(i in 1:cq){
   if(!i %in% unique(pairs[1,close]) ){
    toPrint<-rbind(toPrint,matrix(c(chr,q[i,],rep(NA,ehc),"FN","" ),ncol=length(header)))
   }
  }
  t(toPrint)
 }else if(cq>0){
  #False Negative
  matrix(apply(q,1,function(x) matrix(c(chr,x,rep(NA,ehc),"FN","" ),ncol=length(header))),ncol=length(header),byrow=TRUE)
 }else if(ce>0){
  #False Positive
  t(matrix(apply(e,1,function(x) matrix(c(chr,rep(NA,qhc),x,"FP","" ),ncol=length(header))),ncol=length(header),byrow=TRUE))
 }else{
  #True Negative
  t(matrix(c(chr,rep(NA,qhc+ehc),"TN","" ),ncol=length(header)))
 }
})


if(is.list(toPrint)){
 toPrint<-do.call(rbind,sapply(toPrint,t))
}else{
 toPrint<-matrix(toPrint,ncol=length(header),byrow=TRUE)
}


colnames(toPrint)<-header

write.table(toPrint,file=paste(outprefix,".qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F,na="")

#if(length(estimates[1,])>0){
# colnames(estimates)<-c("chr","start","stop","freqDiff","roughEst")
# write.table(estimates,file=paste(outprefix,".qtlEstimates.csv",sep=""),row.names=F,sep="\t",quote=F)
#}

q("no",0,F)


























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























