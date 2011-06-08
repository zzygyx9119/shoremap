#returns the reached p-value for the interval along with the positions (start, stop, p.value)
#if the third value is larger than the confidence level, it is not a valid region.
#616 not a valid peak
#617 the examined window was too small compared to the internal windowsize (minWindow)

#chromosome,positions, background_count,forground_count and error_count are vectors of the same length
require(bbmle)
require(EMT)

ShoreMap.confint <- function(chromosome,positions,background_count,foreground_count, error_count,foreground_frequency=1,level=0.99,recurse=FALSE,forceInclude=FALSE,allowAdjustment=0.05,filterOutliers=200000,filterPValue=0.05) {
# print(sapply(ls(all.names=TRUE),function(x) eval(parse(text=paste("length(",x,")",sep="")))))
 foreground_frequency<-as.numeric(foreground_frequency)
 print(paste("analysing chr ",chromosome[1],", with ",length(chromosome)," markers for equality to ",foreground_frequency,"(",typeof(foreground_frequency),")",sep=""))
 internalData<- cbind(chromosome,positions,foreground_count,background_count,error_count)
 internalData<- internalData[rowSums(internalData[,3:5])>0,]

 assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
# assign("savedCalc_shoremapmle",0,".GlobalEnv") 
 #apply filtering here:
 filtered=c();
 if(filterOutliers>0){ #condition
  #filterOutliers is the windowsize to use
  f<-filterSampling(internalData,as.numeric(filterOutliers),as.numeric(filterPValue),FALSE)
  print(paste("Removed: ",sum(!f)," markers as outliers"))
  filtered<-internalData[!f,2]
  internalData<- internalData[f,]
 }
 assign("dataset_shoremapmle",internalData,".GlobalEnv")
 freqs<-apply(internalData,1,function(x) x[3]/sum(x[3:5]))
 assign("i_shoremapmle",0,".GlobalEnv")
 minWindow<-max(10,floor(length(internalData[,2])/1000))
 bestsize<- ceiling((max(table(sapply(2:length(freqs),function(x) if(freqs[x]==freqs[x-1]){i_shoremapmle}else{assign("i_shoremapmle",i_shoremapmle+1,".GlobalEnv");i_shoremapmle})))+1)/5)*5
 bestsize<-max(bestsize,minWindow)
 print(paste("bestsize:",bestsize))
 if(bestsize<length(dataset_shoremapmle[,1])){
  ps_global<-sapply(1:(length(dataset_shoremapmle[,1])-bestsize+1),function(x) {cur_freqs<-freqs[x:(x+bestsize-1)]; p<-tryCatch(t.test(cur_freqs,mu=foreground_frequency)$p.value,error=function(err) { -616});ifelse(is.na(p),-616,p)})
  print(paste("finding peaks.. best p-value:",max(ps_global)))
  res<- identify_peaks(1,length(internalData[,2]),foreground_frequency,level,minWindow,ps_global,bestsize,recurse,forceInclude,allowAdjustment)
# rm(dataset_shoremapmle)
# rm(storage_shoremapmle)
  ci<-apply(res,1,function(x) t(c(start=ifelse(x[3]<0,internalData[x[1],2],0), stop=ifelse(x[3]<0,internalData[x[1]+x[2]-1,2],0),p.value=ifelse(x[3]<0,-1*(x[3]+x[2]),x[3]))))
  print("found intervals:")
  apply(ci,2,function(x) print(paste(x[1],"-",x[2])))
  rm(storage_shoremapmle,dataset_shoremapmle,i_shoremapmle)
  list(confidenceInterval=ci,excluded=filtered)
 }else{
  #too few markers
  rm(storage_shoremapmle,dataset_shoremapmle,i_shoremapmle)
  list(confidenceInterval=c(0,0,919),excluded=filtered)
 }
}

filterSampling <- function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE){
 assign("dataset_shoremapmle",internalData,".GlobalEnv")
 fs_ret<- c()
 fs_windows<-floor(internalData[,2]/fs_windowsize)
 fs_allIndices<-1:length(fs_windows)
 print(paste("Analysing ", length(unique(fs_windows)), " ", fs_windowsize, " bp windows for outliers",sep="" ))
 for(fs_window in unique(fs_windows)){
  fs_curIndices<-fs_allIndices[fs_windows==fs_window]
  fs_startPos<-min(fs_curIndices)
  fs_size<-length(fs_curIndices)
  if(fs_size>3){
   fs_p.win<- samplefreqs(fs_startPos,fs_size)
   fs_p.res<-c()
   if(fs_exact){
    sink("/dev/null");
    fs_p.res<-apply(internalData[fs_curIndices,],1,function(x) multinomial.test(x[3:5],prob=fs_p.win)$p.value);
    sink();
   }else{
    fs_p1<-fs_p.win[3]
    fs_p2<-fs_p.win[1]/sum(fs_p.win[1:2])
    fs_p.res<-apply(internalData[fs_curIndices,3:5],1,function(x) pbinom(x[3],size=sum(x),prob=fs_p1)*pbinom(x[1],size=sum(x[1:2]),prob=fs_p2));
   }
   fs_ret<-c(fs_ret,fs_p.res)
  }else{
   fs_ret<-c(fs_ret,rep(1,fs_size))
  }
 }
 p.adjust(fs_ret,method="holm")>=fs_limit
}

filterFreqs <- function(freqs,flanksize, limit=1e-15){
 #calculate adjusted p-values
 outlier_p<-p.adjust(sapply((flanksize+1):(length(freqs)-flanksize),function(x) t.test(c(freqs[(x-flanksize):(x-1)],freqs[(x+1):(x+flanksize)]),mu=freqs[x])$p.value))
 #flank and return
 c(rep(TRUE,flanksize),outlier_p>=limit,rep(TRUE,flanksize))
}

identify_peaks <- function(indexL,indexH,frequency,level,minWindow,ps_global,bestsize,recurse,forceInclude=FALSE,allowAdjustment=0.05){
 if(indexH-indexL>max(minWindow,bestsize)){
#  print(paste(indexL,indexH,bestsize,sep=" ### "))
  cur_indices<-indexL:(indexH-bestsize+1)
  ps<-ps_global[cur_indices]
  if(max(ps)>0.0000000001){
   #try to find peaks
   starts<- cur_indices[ps==max(ps)]
   while(length(starts)>0){
    beststart<- starts[ceiling(length(starts)/2)]
    #see if the frequency needs adjustment
    newFreq<-multll(beststart,bestsize)@coef[1]
    if(abs(frequency-newFreq)<allowAdjustment){
     print(paste("The goal frequency was adjusted from ",frequency," to ",newFreq,", which is the estimated frequency in the preliminary interval",sep=""))
     frequency<-newFreq
    }
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
      resL<- identify_peaks(indexL, res[1,1], frequency, level, minWindow, ps_global, bestsize, recurse, forceInclude,0.0)
      resH<- identify_peaks(res[1,1]+res[1,2], indexH, frequency, level, minWindow, ps_global, bestsize, recurse, forceInclude,0.0)
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

loglikelihood_mult <- function(llm_P1=0.5,llm_err=0.01,llm_index=0,llm_size=0){
 #the loglikelihood function. Returns 110000 for unvalid p-values
# if(P1<0 || P1>1 || err<1e-2 || err>1) {
 if(is.na(llm_P1) || is.na(llm_err)){
  220000
 } else if(llm_P1<0 || llm_P1>1 || llm_err<0 || llm_err>1) {
  110000
 }else{
  llm_P1<-as.numeric(llm_P1)
  llm_err<-as.numeric(llm_err)
  #print(paste(P1,err,sep="##"))
  llm_p1<- llm_P1*(1-4*llm_err/3)+llm_err/3
  llm_pe<- 2*llm_err/3
  llm_p2<- 1-llm_p1-llm_pe
  llm_p.all <- c(llm_p1,llm_p2,llm_pe)
  -sum(apply(dataset_shoremapmle[llm_index:(llm_index+llm_size-1),],1,function(llm_x){dmultinom(x=c(llm_x[3],llm_x[4],llm_x[5]),prob=llm_p.all,log=TRUE)}))
 }
}

samplefreqs <- function(sf_startPos,sf_size) {
 sf_esti<-multll(sf_startPos,sf_size)
 sf_P1<-sf_esti@coef[1]
 sf_err<-sf_esti@coef[2]
 sf_p1<- sf_P1*(1-4*sf_err/3)+sf_err/3
 sf_pe<- 2*sf_err/3
 sf_p2<- 1-sf_p1-sf_pe
 c(sf_p1,sf_p2,sf_pe)
}

multll<- function(ml_x,ml_size=10) {
 #a wrapper for the MLE test
 mle2(loglikelihood_mult,method="Nelder-Mead",start=list(llm_P1=0.5,llm_err=0.0001),fixed=list(llm_index=ml_x,llm_size=ml_size))
}

restrictedModel <- function(P1,x,size) {
 mle2(loglikelihood_mult,optimizer="optimize",start=list(llm_err=0.0001),fixed=list(llm_P1=P1,llm_index=x,llm_size=size),lower=0, upper=1)
}

maxConf<-function(x,level=0.95,freq=0,indexL=0,indexH=Inf,minWindow=10,include=-1){
 #function to minimize for the optimization of the interval
 start<-floor(x[1])
 size<-floor(x[2])
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
   assign("storage_shoremapmle",rbind(storage_shoremapmle,c(start,size,res)),".GlobalEnv")
   res
  }
 }
}

extend <- function(beststart,bestsize=10,level=0.99,freq=1,indexL=0,indexH=Inf,minWindow=10,forceInclude=TRUE){
 #given a window it extends this as far as possible to the left and right without exceeding the confidence level
 bestvalue<-Inf
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