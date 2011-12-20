#lib v01
#v02 - included the absolute height as cutoff for peak prediction
#v02 - fixed a bug in the frequency gradient calculations

library(bbmle)

## interval functions

estimateFreq_one<-function(interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 f<-colSums(shoremap_qtlData[interval,cols])
 f[1]/sum(f)
}


ll_one_fun<-function(f,interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 nrOfPlants<-500
 fc<-min(max(floor(f*nrOfPlants),1),499)/nrOfPlants #fulhack
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==min(interval) & shoremap_qtlmem[,3]==max(interval) & shoremap_qtlmem[,4]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==min(interval) & shoremap_qtlmem[,3]==max(interval)& shoremap_qtlmem[,4]==memory,5]
 }else{
  l<--sum(apply(shoremap_qtlData[interval,cols],1,function(x){dbinom(x[1],size=sum(x),prob=fc,log=TRUE)}))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,min(interval),max(interval),memory,l)),".GlobalEnv")
  l
 }
}

ll_rest_two<-function(f1,f2,minDiff){
 if(f1<0||f1>1){
  110000+abs(f1-0.5)
 }else if(f2<0 || f2>1){
  120000+abs(f2-0.5)
 }else if(abs(abs(f1-f2)-minDiff)>1e-3){
  130000+abs(f1-f2)
 }else{
  ll_one_fun(f1,shoremap_interval,1)+ll_one_fun(f2,shoremap_interval,2)
 }
}


maxConf<-function(x,minDiff,level,validMin=1,validMax=Inf,maxStart=Inf,minStop=-Inf,include=-1){
# print(paste(minDiff,level,sep=" ## "))
 start=floor(x[1])
 stop=floor(x[2])
 validMin<-max(1,validMin)
 validMax<-min(length(shoremap_qtlData[,1]),validMax)
# print(paste(start,stop,sep=" ## "))
 if(start>stop){
  1000000000+start-stop
 }else if(start<validMin){
  1100000000+validMin-start
 }else if(stop>validMax){
  1200000000+stop-validMax
 }else if(start>maxStart){
  1300000000+start-maxStart
 }else if(stop<minStop){
  1400000000+minStop-stop
 }else if(stop-start<3){
  1500000000+start-stop
 }else if(include>0 && !(include>=start && include<=stop)){
  1600000000+start-stop
 }else{
  interval<-start:stop
  f1<-estimateFreq_one(interval,1)
  f2<-estimateFreq_one(interval,2)
  diff<-abs(f1-f2)
  if(abs(minDiff-diff)<1e-3){
   -(shoremap_qtlData[stop,2]-shoremap_qtlData[start,2])
  }else{
   assign("shoremap_interval",interval,".GlobalEnv")
   f1e<-f1
   f2e<-f2
   if(f1e<f2e){
    f1e<-max(0,f1e-(minDiff-diff)/2)
    f2e<-f1e+minDiff
   }else{
    f1e<-min(f1e+(minDiff-diff)/2,1)
    f2e<-f1e-minDiff
   }
   restricted<-mle2(ll_rest_two,start=list(f1=f1e,f2=f2e),fixed=list(minDiff=minDiff))
   full<-ll_one_fun(f1,interval,1)+ll_one_fun(f2,interval,2)
   p<-pchisq(-2*(full-restricted@min),1)
#   print(paste(p,start,stop,sep=" ## "))
   if(p<=level){
    ll<--(shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]+log(1-p))
#    print(paste(p,start,stop,shoremap_qtlData[start,2],shoremap_qtlData[stop,2],ll,sep=" ## "))
    ll
   }else{
    shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]
   }
  }
 }
}


## model functions

singleQtl<-function(q){
# x<-a.t[,1]-q[1]
# tryCatch({
  q[2]*exp(-4/3*(a.t[,1]-q[1])/ifelse(a.t[,1]<q[1],-q[3],q[4]))
# },error=function(err) {
#  cat("err:", conditionMessage(err), "\n")
#  print(q)
#  1
# })
}

freq<-function(x){
 if(length(x)==1){
  #no qtl
  rep(x[1],length(a.t[,1]))
 }else{
  baseline<-x[1]
  qtlest<-matrix(x[-1],byrow=TRUE,ncol=4)
  f<-baseline+rowSums(apply(qtlest,1,singleQtl)) #to handle false frequencies, not caught by the parameter constraints
  f[f>0.99]<-0.99
  f[f<0.01]<-0.01
  f
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
 plot(a.t[,c(1,3)],type="p",pch=16,ylim=c(0,1),col="darkgray",xlab="position",ylab="frequency")
 points(a.t[,1:2],pch=16,col="lightgray")
 if(nrow(qtls)>0){
  apply(matrix(unique(qtls[,2]),ncol=1),1,function(x) abline(v=x))
 }
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
  pm<-which(abs(transitions)==2)
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
 }else if(any(df[c(seq(4,nrow(df),4),seq(5,nrow(df),4)),]<1e7-1)){# optim is c code machine tolerance
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
 }else if(df$x1[1]<df$x2[1]&& (sum(df$x1[seq(3,length(df$x1),4)])>0 ||sum(df$x2[seq(3,length(df$x1),4)])<0)){
  #unbalanced effect and baseline
  if(verbose) print("first")
  13
 }else if(df$x1[1]>df$x2[1]&& (sum(df$x1[seq(3,length(df$x1),4)])<0 ||sum(df$x2[seq(3,length(df$x1),4)])>0)){
  #Failed transition between positive and negative qtls for the second case
  if(verbose) print("second")
  14
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

freqBothGrad_singleQtl<-function(q,d){
 
 s<-sign(a.t[,1]-q[1])
 
 lowIndex<-which(s==-1)
 highIndex<-which(s==1)
 
 low<-length(lowIndex)
 zero<-length(which(s==0))
 high<-length(highIndex)

 (d* singleQtl(q)) %*% cbind(c(rep(-4/3/q[3],low),rep(0,zero),rep(4/3/q[4],high)), rep(1/q[2],length(d)), c(4*(q[1]-a.t[lowIndex,1])/3/q[3]**2,rep(0,zero+high)), c(rep(0,low+zero),4*(a.t[highIndex,1]-q[1])/3/q[4]**2))
}

freqBothGrad<-function(x){
 #analytical gradient of the squared sum
 if(length(x)==2){
  #only for at least 1 QTL
  c(sum(2*(freq(x[1])-a.t[,2])),sum(2*(freq(x[2])-a.t[,3])))
 }else{
  par<-separateParameters(x)
  if(is.na(parameterConstraints(par))){
   NA
  }else{
   d1<-2*(freq(par$x1)-a.t[,2])
   d2<-2*(freq(par$x2)-a.t[,3])
   c(sum(d1), sum(d2), apply(rbind(apply(matrix(par$x1[-1],byrow=TRUE,ncol=4), 1, freqBothGrad_singleQtl,d=d1), apply(matrix(par$x2[-1],byrow=TRUE,ncol=4), 1, freqBothGrad_singleQtl,d=d2)),2,function(col) c(col[1]+col[5],col[2:4],col[6:8])))
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
  -sum(apply(cbind(data[,2],rowSums(data[,2:3]),freq(df$x1)),1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))-sum(apply(cbind(data[,4],rowSums(data[,4:5]),freq(df$x2)),1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))
 }
}

getNewQtl<-function(x,direction=c(-1,1)){
 chrsize<-max(a.t[,1])+1
 posToUse<-0
 effect<-c(0,0)
 if(length(x)==2){
  posToUse<-round(chrsize/2)
  effect<-sign(x-0.5)*1e-2
 }else{
  df<-separateParameters(x)
  #find section with most deviation per marker
  pos<-matrix(sort(df[seq(2,nrow(df),4),1]),ncol=1)
  windows<-rowSums(apply(pos,1,function(p) ifelse(a.t[,1]<p,0,1)))
  d1<-freq(df$x1)-a.t[,2]
  d2<-freq(df$x2)-a.t[,3]
  winToUse<-which.max(tapply(d1**2+(d2)**2,windows,mean))
  posToUse<-round((diff(c(0,pos,chrsize))/2+c(0,pos))[winToUse])
  effect<-sign(c(tapply(d1,windows,mean)[winToUse],tapply(d2,windows,mean)[winToUse]))*1e-2
 }
 #set direction of effects
 effect<-direction*3e-2
 par<-c(posToUse,effect[1],1e7,1e7,effect[2],1e7,1e7)
 if(is.na(freqBoth(c(x,par)))){
  c(posToUse,0,1e7,1e7,0,1e7,1e7)
 }else{
  par
 }
}


calcFocusScale<-function(x){
 len<-length(x)
 if(len>2){
  c(1,1, rep(c(max(a.t[,1])+1, rep(c(2, rep(max(x[c(seq(5,len,7), seq(6,len,7), seq(8,len,7), seq(9,len,7))]),2)),2)), length(x)/7))
 }else{
  c(1,1)
 }
}


calcbic<-function(x){
 2*llBoth(x)+length(x)*log(nrow(data))
}

findBest<-function(par.start,tolerance=1e4,focus=TRUE,verbose=FALSE){
 if(verbose){
  estPlot(par.start)
 }
 scale<-par.start
 if(focus){
  scale<-calcFocusScale(par.start)
 }
 opt<-optim(fn=freqBoth, par=par.start, method="Nelder-Mead", control=list(parscale=scale, maxit=500))
 if(verbose){
  estPlot(opt$par)
  print(opt$val)
 }
 val<-Inf
 while (opt$val<val-tolerance && opt$convergence==1) {
  val<-opt$val
  scale<-opt$par
  if(focus){
   scale<-calcFocusScale(opt$par)
  }
  opt<-optim(fn=freqBoth,gr=freqBothGrad, par=opt$par, method="BFGS", control=list(parscale=scale, maxit=500))
  scale<-opt$par
  if(focus){
   scale<-calcFocusScale(opt$par)
  }
  opt<-optim(fn=freqBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=scale, maxit=500))
  if(verbose){
   estPlot(opt$par)
   print(opt$val)
  }
 }
 opt
}

refine<-function(par.start,tolerance=1e4,focus=TRUE,verbose=FALSE){
 if(verbose){
  print("Refining...")
 }
 scale<-par.start
 if(focus){
  scale<-calcFocusScale(par.start)
 }
 opt<-optim(llBoth,par=par.start, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
 if(verbose){
  estPlot(opt$par)
  print(opt$val)
 }
 val<-Inf
 while (opt$val<val-tolerance && opt$convergence==1) {
 #while (opt$convergence==1) {
  val<-opt$val
  scale<-opt$par
  if(focus){
   scale<-calcFocusScale(opt$par)
  }
  opt<-optim(llBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=scale, maxit=5000))
  if(verbose){
   estPlot(opt$par)
   print(opt$val)
  }
 }
 opt
}


