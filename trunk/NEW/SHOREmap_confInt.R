#returns the reached p-value for the interval along with the positions (start, stop, p.value)
#if the third value is larger than the confidence level, it is not a valid region.
#616 not a valid peak
#617 the examined window was too small compared to the internal windowsize (minWindow)

#chromosome,positions, background_count,forground_count and error_count are vectors of the same length
require(bbmle)

ShoreMap.confint <- function(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,recurse=FALSE,forceInclude=FALSE){
 foreground_frequency<-as.numeric(foreground_frequency)
 print(paste("analysing chr ",chromosome[1],", with ",length(chromosome)," markers for equality to ",foreground_frequency,"(",typeof(foreground_frequency),")",sep=""))
 internalData<- cbind(chromosome,positions,foreground_count,background_count,error_count)
 internalData<- internalData[rowSums(internalData[,3:5])>0,]
 #perhaps not the nicest solution
 assign("dataset_shoremapmle",internalData,".GlobalEnv")
 assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
# assign("savedCalc_shoremapmle",0,".GlobalEnv")
 minWindow=max(10,floor(length(internalData[,2])/1000))
 assign("i_shoremapmle",0,".GlobalEnv")
 freqs<-apply(dataset_shoremapmle,1,function(x) x[3]/sum(x[3:5]))
 bestsize<- ceiling((max(table(sapply(2:length(freqs),function(x) if(freqs[x]==freqs[x-1]){i_shoremapmle}else{assign("i_shoremapmle",i_shoremapmle+1,".GlobalEnv");i_shoremapmle})))+1)/5)*5
 bestsize<-max(bestsize,minWindow)
 print(paste("bestsize:",bestsize))
 if(bestsize<length(dataset_shoremapmle[,1])){
  ps_global<-sapply(1:(length(dataset_shoremapmle[,1])-bestsize+1),function(x) {cur_freqs<-freqs[x:(x+bestsize-1)]; p<-tryCatch(t.test(cur_freqs,mu=foreground_frequency)$p.value,error=function(err) { -616});ifelse(is.na(p),-616,p)})
  print(paste("finding peaks.. best p-value:",max(ps_global)))
  res<- identify_peaks(1,length(internalData[,2]),foreground_frequency,level,minWindow,ps_global,bestsize,recurse,forceInclude)
# rm(dataset_shoremapmle)
# rm(storage_shoremapmle)
  ci<-apply(res,1,function(x) t(c(start=ifelse(x[3]<0,internalData[x[1],2],0), stop=ifelse(x[3]<0,internalData[x[1]+x[2]-1,2],0),p.value=ifelse(x[3]<0,-1*(x[3]+x[2]),x[3]))))
  print("found intervals:")
  apply(ci,2,function(x) print(paste(x[1],"-",x[2])))
  ci
 }else{
  #too few markers
  c(0,0,919)
 }
}

identify_peaks <- function(indexL,indexH,frequency,level,minWindow,ps_global,bestsize,recurse,forceInclude=FALSE){
 if(indexH-indexL>max(minWindow,bestsize)){
#  print(paste(indexL,indexH,bestsize,sep=" ### "))
  cur_indices<-indexL:(indexH-bestsize+1)
  ps<-ps_global[cur_indices]
  if(max(ps)>0.0000000001){
   #try to find peaks
   starts<- cur_indices[ps==max(ps)]
   while(length(starts)>0){
    beststart<- starts[ceiling(length(starts)/2)]
    res<-extend(beststart,bestsize,level,frequency,indexL,indexH,minWindow, forceInclude)
    if(res[1,3]>0){
     if(length(starts)==1){
      return(t(as.matrix(c(1,1,616))))
      break
     }else{
      starts<-starts[1:length(starts)!=ceiling(length(starts)/2)]
     }
    }else if(recurse){
     if(res[1,3]<0){
      resL<- identify_peaks(indexL, res[1,1], frequency, level, minWindow, ps_global, bestsize, recurse, forceInclude)
      resH<- identify_peaks(res[1,1]+res[1,2], indexH, frequency, level, minWindow, ps_global, bestsize, recurse, forceInclude)
      if(resL[1,3]<0){
       if(resH[1,3]<0){
        #both good
        return(rbind(resL,res,resH))
        break
       }else{
        #low good
        return(rbind(resL,res))
        break
       }
      }else if(resH[1,3]<0){
       #high good
       return(rbind(res,resH))
        break
      }else{
       return(res)
       break
      }
     }else{
      return(res)
      break
     }
    }else{
     return(res)
     break
    }
   }#end while loop
  }else{
   #No peak
   return(t(as.matrix(c(1,1,616))))
  }
 }else{
  #Too small window
  return(t(as.matrix(c(1,1,617))))
 }
}

loglikelihood_mult <- function(P1=0.5,err=0.01,index=0,size=0){
 #the loglikelihood function. Returns 110000 for unvalid p-values
# if(P1<0 || P1>1 || err<1e-2 || err>1) {
 if(is.na(P1) || is.na(err)){
  220000
 } else if(P1<0 || P1>1 || err<0 || err>1) {
  110000
 }else{
  P1=as.numeric(P1)
  err=as.numeric(err)
  #print(paste(P1,err,sep="##"))
  p1<- P1*(1-4*err/3)+err/3
  pe<- 2*err/3
  p2<- 1-p1-pe
  p.all <- c(p1,p2,pe)
  -sum(apply(dataset_shoremapmle[index:(index+size-1),],1,function(x){dmultinom(x=c(x[3],x[4],x[5]),prob=p.all,log=TRUE)}))
 }
}

multll<- function(x,size=10) {
# if(x%%100==0){
#  print(x)
# }
 #a wrapper for the MLE test
 mle2(loglikelihood_mult,method="Nelder-Mead",start=list(P1=0.5,err=0.0001),fixed=list(index=x,size=size))
}

restrictedModel <- function(P1,x,size) {
 mle2(loglikelihood_mult,optimizer="optimize",start=list(err=0.0001),fixed=list(P1=P1,index=x,size=size),lower=0, upper=1)
}

maxConf<-function(x,level=0.95,freq=0,indexL=0,indexH=Inf,minWindow=10,include=-1){
 #function to minimize for the optimization of the interval
 start<-floor(x[1])
 size<-floor(x[2])
# print(paste("start:",start,", size:",size))
 indexL<- max(1,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 if(size<minWindow){
  110000-size+minWindow
 }else if(start<indexL){
  110000+indexL-start
 }else if(start+size-1>indexH){
  110000-indexH-1+start+size
 }else if(include!=-1 && start>include){
  #given point not included
  110000+start-include
 }else if(include!=-1 && start+size-1<include){
  #given point not included
  110000+include-start-size+1
 }else{
  #check storage
  if(sum(storage_shoremapmle[,1]==start&storage_shoremapmle[,2]==size)==1){
   res<-storage_shoremapmle[storage_shoremapmle[,1]==start&storage_shoremapmle[,2]==size,3]
#   print(paste(start,size,res,sep=" # "))
#   assign("savedCalc_shoremapmle",savedCalc_shoremapmle+1,".GlobalEnv")
   res
  }else{
   #if not in storage, calculate
   fit<-multll(start,size)
#   fit<- restrictedModel(freq-1/1000,start,size)
   if(fit@min>100000){
    res<-fit@min
   }else{
    restrictedFit<- restrictedModel(freq,start,size)
#    p<- pchisq(-2*(fit@min-restrictedFit@min),1)
    p<- pchisq(-2*(loglikelihood_mult(fit@coef[1],fit@coef[2],start,size)-loglikelihood_mult(freq,restrictedFit@coef[1],start,size)),1)
    if(p<=level){
     res<- -size-p
    }else{
     res<- size+p
    }
   }
#   print(paste(start,size,res,sep="   "))
   assign("storage_shoremapmle",rbind(storage_shoremapmle,c(start,size,res)),".GlobalEnv")
   res
  }
 }
}

extend <- function(beststart,bestsize=10,level=0.99,freq=1,indexL=0,indexH=Inf,minWindow=10,forceInclude=TRUE){
 #given a window it extends this as far as possible to the left and right without exceeding the confidence level
 bestvalue=Inf
 indexL<- max(1,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 inclusion<--1
 if(forceInclude){
  inclusion<-beststart+floor(bestsize/2)
 }
 
 #a first optimization
 nextTest<- optim(fn=maxConf,par=c(beststart,bestsize),control=list(ndeps=c(1,1)),level=level,freq=freq,indexL=max(indexL,beststart-10*minWindow),indexH=min(indexH,beststart+bestsize+10*minWindow),minWindow=minWindow,include=inclusion)



 while(nextTest$value<bestvalue){
  bestvalue<-nextTest$value
  beststart<-floor(nextTest$par[1])
  bestsize<-floor(nextTest$par[2])
  #alternatively right-left extension by minWindow
  i<-1
  lastvalue<-bestvalue
  curvalue<-maxConf(c(beststart-floor(i/2)*minWindow,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion) 
  while(i<1000&&curvalue<lastvalue){
   lastvalue<-curvalue
   i<-i+1
   curvalue<-maxConf(c(beststart-floor(i/2)*minWindow,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
  }
  i<-i-1
  beststart<-beststart-floor(i/2)*minWindow
  bestsize<-bestsize+i*minWindow
  if(i%%2==0){
   #if the extension breaks down on the right, continue on to expand left
   i<-1
   lastvalue<-maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   curvalue<-maxConf(c(beststart-i*minWindow,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart-i*minWindow,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   }
   i<-i-1
   beststart<-beststart-i*minWindow
   bestsize<-bestsize+i*minWindow
  }else{
   #if the extension breaks down on the left, continue on to expand right
   i<-1
   lastvalue<-maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   curvalue<-maxConf(c(beststart,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart,bestsize+i*minWindow),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   }
   i<-i-1
   beststart<-beststart
   bestsize<-bestsize+i*minWindow
  }
  #alternatively right-left extension
  i<-1
  lastvalue<- maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
  curvalue<-maxConf(c(beststart-floor(i/2),bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
  
  while(i<1000&&curvalue<lastvalue){
   lastvalue<-curvalue
   i<-i+1
   curvalue<-maxConf(c(beststart-floor(i/2),bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
  }
  i<-i-1
  beststart<-beststart-floor(i/2)
  bestsize<-bestsize+i
  if(i%%2==0){
   #if the extension breaks down on the right, continue on to expand left
   i<-1
   lastvalue<-maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   curvalue<-maxConf(c(beststart-i,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart-i,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   }
   i<-i-1
   beststart<-beststart-i
   bestsize<-bestsize+i
  }else{
   #if the extension breaks down on the left, continue on to expand right
   i<-1
   lastvalue<-maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   curvalue<-maxConf(c(beststart,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   while(i<1000&&curvalue<lastvalue){
    lastvalue<-curvalue
    i<-i+1
    curvalue<-maxConf(c(beststart,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   }
   i<-i-1
   beststart<-beststart
   bestsize<-bestsize+i
  }
  #optimization again before reitteration
  nextTest<- optim(fn=maxConf,method="Nelder-Mead",par=c(beststart,bestsize),control=list(ndeps=c(1,1),maxit=100),level=level,freq=freq,indexL=max(indexL,beststart-10*minWindow),indexH=min(indexH,beststart+bestsize+10*minWindow),minWindow=minWindow,include=inclusion)
 }
 t(as.matrix(c(beststart,bestsize,bestvalue)))
}
