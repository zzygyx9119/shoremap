#chromosome,positions, background_count,foreground_count and error_count are vectors of the same length
ShoreMap.confint <- function(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,minWindow=10,recurse=TRUE){
 internalData<- cbind(chromosome,positions,foreground_count,background_count,error_count)
 internalData<- internalData[rowSums(internalData[,3:5])>0,]
 #perhaps not the nicest solution
 assign("dataset_shoremapmle",internalData,".GlobalEnv")
 res<- identify_peaks(0,length(internalData[,2]),foreground_frequency,level,minWindow,recurse)
 apply(res,1,function(x) t(c(start=internalData[x[1],2], stop=internalData[x[1]+x[2]-1,2],p.value=-1*(x[3]+x[2]))))
}

identify_peaks <- function(indexL,indexH,frequency,level,minWindow,recurse){
 beststart=floor(optimize(maxConfPos,lower=indexL,upper=indexH,freq=frequency,neighbors=floor(minWindow/4),size=minWindow,minWindow=minWindow)$minimum[1])
# might avoid a check by checking the first region
# maxConf(c(beststart,10),level=level,freq=frequency,indexL=indexL,indexH=indexH,minWindow=minWindow)
 res<-extend(beststart,minWindow,level,frequency,indexL,indexH,minWindow)
 if(res[1,3]<0){
  #try to find other peaks to the left and right
  if(recurse){
   resL<- identify_peaks(indexL,max(res[1,1],indexL,0),frequency,level,minWindow)
   resH<- identify_peaks(min(res[1,1]+res[1,2],indexH,length(dataset_shoremapmle[,2])),indexH,frequency,level,minWindow)
   if(resL[1,3]<0){
    if(resH[1,3]<0){
     #both good
     rbind(resL,res,resH)
    }else{
     #low good
     rbind(resL,res)
    }
   }else if(resH[1,3]<0){
    #high good
    rbind(res,resH)
   }else{
    res
   }
  }else{
   res
  }
 }else{
  res
 }
}

loglikelihood_mult <- function(P1=0.5,err=0.01,index=0,size=0){
 #the loglikelihood function. Returns 110000 for unvalid p-values
 if(P1<0 || P1>1 || err<0 || err>1) {
  110000
 }else{
  p1<- P1*(1-4*err/3)+err/3
  pe<- 2*err/3
  p2<- 1-p1-pe
  p.all <- c(p1,p2,pe)
  -sum(apply(dataset_shoremapmle[index:(index+size-1),],1,function(x){dmultinom(x=c(x[3],x[4],x[5]),prob=p.all,log=TRUE)}))
 }
}

multll<- function(x,size=10) {
 #a wrapper for the MLE test
 mle2(loglikelihood_mult,method="Nelder-Mead",start=list(P1=0.5,err=0.0001),fixed=list(index=x,size=size))
}

maxConf<-function(x,level=0.95,freq=0,indexL=0,indexH=Inf,minWindow=10){
 #function to minimize for the optimization of the interval
 start<-floor(x[1])
 size<-floor(x[2])
 indexL<- max(0,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 if(size<minWindow){
  110000-size+minWindow
 }else if(start<indexL){
  110000+indexL-start
 }else if(start+size-1>indexH){
  110000+indexH+1-start-size
 }else{
  fit<-multll(start,size)
  if(fit@min>100000){
   res<-fit@min
  }else{
   p<- pchisq(-2*(loglikelihood_mult(fit@coef[1],fit@coef[2],start,size)-loglikelihood_mult(freq,fit@coef[2],start,size)),1)
   if(p<=level){
    res<- -size-p
   }else{
    res<- size+p
   }
  }
#  print(paste(start,size,res,sep="   "))
  res
 }
}


maxConfPos<-function(x,level=0.95,freq=0,size=10,neighbors=2,indexL=0,indexH=Inf,minWindow=10){
 #function to optimize in order to find an initial start point
 start<-floor(x[1])
 indexL<- max(0,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 if(start<indexL+neighbors){
  110000-indexL-neighbors+start
 }else if(size<minWindow){
  110000-size+minWindow
 }else if(start+size-1>indexH-neighbors){
  110000+start+size-1-indexH+neighbors
 }else{
  res<-sum(sapply((start-2):(start+2),function(y) abs(multll(y,size)@coef[1]-freq)))
  res
 }
}




extend <- function(beststart,bestsize=10,level=0.95,freq=1,indexL=0,indexH=Inf,minWindow){
 #given a window it extends this as far as possible to the left and right without exceeding the confidence level
 bestvalue=Inf
 indexL<- max(0,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 #a first optimization
 nextTest<- optim(fn=maxConf,par=c(beststart,bestsize),control=list(ndeps=c(1,1)),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)



 while(nextTest$value<bestvalue){
  bestvalue<-nextTest$value
  beststart<-floor(nextTest$par[1])
  bestsize<-floor(nextTest$par[2])
  i<-1
  lastvalue<-bestvalue
  curvalue<-maxConf(c(beststart-floor(i/2),bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
  #alternatively right-left extension
  while(i<1000&&curvalue<lastvalue){
   lastvalue<-curvalue
   i<-i+1
   curvalue<-maxConf(c(beststart-floor(i/2),bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
  }
  i<-i-1
  beststart<-beststart-floor(i/2)
  bestsize<-bestsize+i
  if(i%%2==0){
   #if the extension breaks down on the right, continue on to expand left
   i<-1
   lastvalue<-bestvalue
   curvalue<-maxConf(c(beststart-i,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart-i,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
   }
   i<-i-1
   beststart<-beststart-i
   bestsize<-bestsize+i
  }else{
   #if the extension breaks down on the left, continue on to expand right
   i<-1
   lastvalue<-bestvalue
   curvalue<-maxConf(c(beststart,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
   }
   i<-i-1
   beststart<-beststart
   bestsize<-bestsize+i
  }
  #optimization again before reitteration
  nextTest<- optim(fn=maxConf,par=c(beststart,bestsize),control=list(ndeps=c(1,1),maxit=100),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow)
 }
 t(as.matrix(c(beststart,bestsize,bestvalue)))
}
