#v01 - combining qtlseq v13 and testModel v01. modified the interval-finding to find an interval that matches the frequency at the point estimator within 1e-3. let's see if that works
#v02 - small bugfixes
library(bbmle)
source("/projects/dep_coupland/grp_nordstrom/projects/Lotus/Simulation/PopSimulatorRIL/Rcode/SHOREmap_qtlseq_withModel_lib.R")

#Run.....

args<-commandArgs(trailingOnly=TRUE)

#command line arguments
high_file<-args[1]
low_file<-args[2]

if(length(args)==4){
 qtl_file<-args[3]
 
 outprefix<-args[4]
}else if(length(args)==3){
 outprefix<-args[3]
 qtl_file=""
}else{
 print("Usage:\nR --args seqfile1 seqfile2 outprefix < SHOREmap_qtlseq.R\nor\nR --args seqfile1 seqfile2 qtlfile outprefix < SHOREmap_qtlseq.R")
 q("no",616,FALSE)
}


#Parameters
tolerance=2e6
winSize<-1000000
winStep<-10000
p.valueCutoff<-0.2
minMarkersInWindow<-10
minCoverage <- 0
maxCoverage <- 3000000
minFrequencyChange<-0.05
minAbsolutePeak<-0.05
bootstrap=1000

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



qtls<-data.frame(pos=-1,chr=-1)
if(qtl_file==""){
 qtls<-data.frame(pos=-1,chr=-1)
}else{
 qtls<-read.table(paste(qtl_file,".qtl.txt",sep=""))
 ls.t<-read.table(paste(qtl_file,".low.AF.txt",sep=""))
 hs.t<-read.table(paste(qtl_file,".high.AF.txt",sep=""))
 a.t<-cbind(ls.t[,1],ls.t[,2],(hs.t[,3]-ls.t[,3])/100)
 marker<-qtls[,2]+1
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
 rank<-sort(sort(abs(qtls[,5]),index.return=TRUE,decreasing=TRUE)$ix,index.return=TRUE)$ix
 lg<-sapply(qtls[,3],function(x) ifelse(sum(qtls[,3]==x)>1,x,0))
 lgType<-sapply(lg,function(x) ifelse(x==0,0,{
  tt<-table(sign(qtls[lg==x,5]))
  ifelse(length(tt)==1,sign(sum(qtls[lg==x,5])),ifelse(tt[1]<tt[2],-tt[1]/tt[2],ifelse(tt[1]==tt[2],0,tt[2]/tt[1])))
 }))
 qtls<-data.frame(id=qtls[,1],pos=qtls[,4],chr=qtls[,3],effect=qtls[,5],rank=rank,lg=lg,lgType=lgType,trueFreqDiff=freqDiff,trueStart=trueStart,trueEnd=trueEnd)
 
}


##merge
hs.freq<-hs[,3]/rowSums(hs[,3:5])
ls.freq<-ls[,3]/rowSums(ls[,3:5])

data_all<- cbind(hs[,1:2],hs.freq,hs[,3:5],ls.freq,ls[,3:5],hs.freq-ls.freq)

#p<-apply(data,1,function(x) fisher.test(matrix(as.numeric(x[c(4,5,8,9)]),nrow=2))$p.value)
p<-rep(0,nrow(data_all))

data_all<-cbind(data_all,p)


estimates<- matrix(ncol=5)[FALSE,]
colnames(estimates)<-c("chr","start","stop","freqDiff","Est")

nrOfChrs<-length(unique(data_all[,1]))
verbose<-FALSE #the implementation needs to be changed
pdf(paste(outprefix,".plots.pdf",sep=""))
par(mfrow=c(ceiling(nrOfChrs/2),2))

colors<-rainbow(50,start=4/6,end=1/6)

for(chr in unique(data_all[,1])){
 print(paste("Chr:",chr))
 #data for interval
 data2<-data_all[data_all[,1]==chr,]
 #data for model
 a.t<-t(sapply(data2[,2],function(x) {index<-which(data2[,2]>=x-winSize/2 & data2[,2]<x+winSize/2);c(x,sum(data2[index,8])/sum(data2[index,8:9]),sum(data2[index,4])/sum(data2[index,4:5]))}))
 data<- data2[,c(2,8,9,4,5)]


# png()
 plot(data2[,2],data2[,11],main=paste("chr",chr),xlab="pos",ylab="frequency difference",ylim=c(-1.2,1),pch=16,cex=0.75,col="lightsteelblue3")
 lines(a.t[,1],a.t[,3]-a.t[,2],col="limegreen")
 #calculate peaks
 if(verbose) X11()
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
 #output parameters to file
 sapply(pars,function(x) write.table( matrix(c(chr,calcbic(x$par),x$par),nrow=1),file=paste(outprefix,"_modelParameters_raw.csv",sep=""), append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE))

 bestPar<-pars[[which.min(sapply(pars,function(x) calcbic(x$par)))]]$par
 if(TRUE){ #add refinement by further taking binomial sampling into account
  bestPar<-refine(bestPar,verbose=verbose)$par
 }
 write.table( matrix(c(chr,calcbic(bestPar),bestPar),nrow=1),file=paste(outprefix,"_modelParameters_refined.csv",sep=""), append=TRUE, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

 if(verbose) dev.off()
 
 if(length(bestPar)>2){
  df<-separateParameters(bestPar)
  assign("shoremap_qtlData",data2,".GlobalEnv")  
  for(peakIndex in seq(2,length(df$x1),4)){
   assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1,-1),nrow=1),".GlobalEnv")
   markerIndex<-which.min(abs(data2[,2]-df$x1[peakIndex]))
   markerDiff<-abs(a.t[markerIndex,2]-a.t[markerIndex,3])
   startIndex<-max(1,markerIndex-5)
   endIndex<-min(nrow(data2),markerIndex+5)



  

    #identify window
   opt<-optim(fn=maxConf,par=c(startIndex,endIndex),minDiff=markerDiff,level=0.99,include=markerIndex)
   bestValue<-Inf
   while(opt$value<bestValue){
    bestValue<-opt$value
    opt<-optim(fn=maxConf,par=c(floor(opt$par[1]),floor(opt$par[2])),minDiff=markerDiff,level=0.99,include=markerIndex)
   }
   rect(data2[opt$par[1],2],0.96,data2[opt$par[2],2],1,col="limegreen",border="limegreen")

   est<-df$x1[peakIndex]
   abline(v=est,col="limegreen")

   estimates<-rbind(estimates,c(chr,data2[opt$par[1],2],data2[opt$par[2],2],markerDiff,est))
  }
 }
 if(sum(qtls$chr==chr)>0){
  apply(qtls[qtls$chr==chr,],1,function(x) {
   e<-min(50,x[4])
   e<-max(-50,e)
   e<-round((e+50)/2)
   abline(v=x[2],col=colors[e])
  })
 }
 #print scale
 chrStart<-min(data2[,2])
 chrEnd<-max(data2[,2])
 chrSize<-chrEnd-chrStart
 scaleStart<-chrStart+chrSize/5
 scaleEnd<-chrEnd-chrSize/5
 scaleSize<-scaleEnd-scaleStart
 dx<-scaleSize/50
 starts<-seq(scaleStart,scaleEnd-1,dx)  
  
 rect(starts, -1, starts+0.9*dx, -1.1, col = colors, border = "grey");
 text(scaleStart+scaleSize/2,-0.95,"effects")
 text(scaleStart+scaleSize/2,-1.15,"0")
 text(starts[1],-1.15,"-50")
 text(max(starts)+0.9*dx,-1.15,"50")

}

dev.off()


#print estimates
header<-c("chr",colnames(qtls)[c(1:2,4:ncol(qtls))],colnames(estimates)[c(2:ncol(estimates))],"judgement","spec")
header[1]<-paste("#",header[1],sep="")
qhc<-ncol(qtls)-1
ehc<-ncol(estimates)-1

toPrint<-sapply(unique(data_all[,1]),function(chr){
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
  t(matrix(apply(q,1,function(x) matrix(c(chr,x,rep(NA,ehc),"FN","" ),ncol=length(header))),ncol=length(header),byrow=TRUE))
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
