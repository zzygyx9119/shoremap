#library v01

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
 plot(a.t[,c(1,3)],type="l",ylim=c(0,1),xlab="position",ylab="frequency")
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
 effect<-sign(direction)*5e-2
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

addQtl<-function(previous,verbose=FALSE){
 par.start<-getNewQtl(previous)
 if(is.na(parameterConstraints(separateParameters(c(previous,par.start))))){
  par.start<-getNewQtl(previous,c(1,-1))
 }
 if(verbose){
  estPlot(c(previous,par.start))
 }
 opt<-optim(freqBoth_addOne, par=par.start, static=previous, method="Nelder-Mead", control=list(parscale=par.start, maxit=5000))
 if(verbose){
  estPlot(c(previous,opt$par))
  print(opt$val)
 }
 val<-Inf
 while (opt$val<val-1e-3 && opt$convergence==1) {
  val<-opt$val
  opt<-optim(freqBoth, par=opt$par, static=previous, method="Nelder-Mead", control=list(parscale=opt$par, maxit=5000))
  if(verbose){
   estPlot(c(previous,opt$par))
   print(opt$val)
  }
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

freqBothGrad_singleQtl<-function(q,d){
 
 count<-table(sign(a.t[,1]-q[1]))
 low<-0
 zero<-0
 high<-0
 if(any(names(count)=="-1")){
  low<-count[names(count)=="-1"]
 }
 if(any(names(count)=="0")){
  zero<-count[names(count)=="0"]
 }
 if(any(names(count)=="1")){
  high<-count[names(count)=="1"]
 }
 (d* singleQtl(q)) %*% cbind(c(rep(-4/3/q[3],low),rep(0,zero),rep(4/3/q[4],high)), rep(1/q[2],length(d)), c(4*(q[1]-a.t[1:low,1])/3/q[3]**2,rep(0,zero+high)), c(rep(0,low+zero),4*(a.t[(length(d)-high+1):length(d),1]-q[1])/3/q[4]**2))
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
#   qtlest1<-matrix(par$x1[-1],byrow=TRUE,ncol=4)
#   qtlest2<-matrix(par$x2[-1],byrow=TRUE,ncol=4)
#   qtlGrad1<-apply(qtlest1,1,freqBothGrad_singleQtl,d=d1)
#   qtlGrad2<-apply(qtlest2,1,freqBothGrad_singleQtl,d=d2)
#   c(sum(d1), sum(d2), apply(rbind(qtlGrad1,qtlGrad2),2,function(col) c(col[1]+col[5],col[2:4],col[6:8])))
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

generateConstraints<-function(x){
 chrsize<-max(a.t[,1])+1
 len<-length(x)
 #b1>0
 ui<-c(1,rep(0,len-1))
 ci<-0
 #b1<1
 ui<-rbind(ui,c(-1,rep(0,len-1)))
 ci<-rbind(ci,-1)
 #b2>0
 ui<-rbind(ui,c(0,1,rep(0,len-2)))
 ci<-rbind(ci,0)
 #b2<1
 ui<-rbind(ui,c(0,-1,rep(0,len-2)))
 ci<-rbind(ci,-1)
 #0<pos<chrsize
 for(i in seq(3,len,7)){
  constr<-rep(0,len)
  constr[i]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,0)
  constr[i]<--1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-chrsize)
 }
 #-1<effect<1 && 0<effect+b<1
 for(i in seq(4,len,7)){
  constr<-rep(0,len)
  constr[i]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)
  constr<-rep(0,len)
  constr[i]<--1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)
  constr<-rep(0,len)
  constr[i]<-1
  constr[1]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,0)
  constr[i]<--1
  constr[1]<--1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)  
 }
 for(i in seq(7,len,7)){
  constr<-rep(0,len)
  constr[i]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)
  constr<-rep(0,len)
  constr[i]<--1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)
  constr<-rep(0,len)
  constr[i]<-1
  constr[2]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,0)
  constr[i]<--1
  constr[2]<--1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,-1)  
 }
 #sigma>5e6
 for(i in c(seq(5,len,7),seq(6,len,7),seq(8,len,7),seq(9,len,7))){
  constr<-rep(0,len)
  constr[i]<-1
  ui<-rbind(ui,constr)
  ci<-rbind(ci,5e6)
 }
 rownames(ui)<-NULL
 rownames(ci)<-NULL
 list(ui=ui,ci=ci)
}

#calcFocusScale<-function(x){
# if(length(x)>2){
#  x*c(1,1,rep(0.5,7),rep(1,length(x)-9))
# }else{
#  x
# }
#}

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
 if(verbose){
  print("Refining...")
 }
 par.start<-opt$par
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


