#returns the reached p-value for the interval along with the positions (start, stop, p.value)
#if the third value is larger than the confidence level, it is not a valid region.
#616 not a valid peak
#617 the examined window was too small compared to the internal windowsize (minWindow)
#919 too few markers after removal of 0-sum-markers, filtering. Too few means that there are less than ten or that next too all are constant

#chromosome,positions, background_count,forground_count and error_count are vectors of the same length
require(bbmle)
require(EMT)

ShoreMap.confint <- function(chromosome,positions, background_count, foreground_count, error_count, foreground_frequency=1, level=c(0.95,0.99,0.999), recurse=FALSE, forceInclude=TRUE, allowAdjustment=0.0, filterOutliers=200000, filterPValue=0.05,winSize=50000,minMarker=0,minCoverage=0) {
# allowAdjustment=0.0
# minMarker=10
# print(sapply(ls(all.names=TRUE),function(x) eval(parse(text=paste("length(",x,")",sep="")))))
 foreground_frequency<-as.numeric(foreground_frequency)
 internalData<- cbind(chromosome,positions,foreground_count,background_count,error_count)
 internalData<- internalData[rowSums(internalData[,3:5])>minCoverage,]
 print(paste("Analysing chr ",chromosome[1],", with ",length(chromosome)," (",length(internalData[,1]),") markers for equality to ",foreground_frequency,"(",typeof(foreground_frequency),")",sep=""))
 assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
# assign("savedCalc_shoremapmle",0,".GlobalEnv") 
 #apply filtering here:
 filtered=c();
 if(filterOutliers>0){ #condition
  #filterOutliers is the windowsize to use
  f<-filterSamplingv3(internalData,as.numeric(filterOutliers),as.numeric(filterPValue),FALSE)
  print(paste("Removed: ",sum(!f)," markers as outliers"))
  filtered<-internalData[!f,2]
  internalData<- internalData[f,]
 }
 assign("dataset_shoremapmle",internalData,".GlobalEnv")
 freqs<-internalData[,3]/rowSums(internalData[,3:5])
 assign("i_shoremapmle",0,".GlobalEnv")
 minWindow<-max(10,floor(length(internalData[,2])/10000))
 bestsize<- ceiling((max(table(sapply(2:length(freqs),function(x) if(freqs[x]==freqs[x-1]){i_shoremapmle}else{assign("i_shoremapmle",i_shoremapmle+1,".GlobalEnv");i_shoremapmle})))+1)/5)*5
 bestsize<-max(bestsize,minWindow)
 print(paste("Bestsize:",bestsize))
 if(bestsize<length(dataset_shoremapmle[,1])){
  avg_pos<-c(sapply(seq(0,winSize-1,10000),function(shift){
   windows<- floor((internalData[,2]+shift)/winSize)
   windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
   tapply(internalData[windowsToUse,2],windows[windowsToUse],mean)
  }),recursive=TRUE)
  avg_freq<-c(sapply(seq(0,winSize-1,10000),function(shift){
   windows<- floor((internalData[,2]+shift)/winSize)
   windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
   tapply(internalData[windowsToUse,3],windows[windowsToUse],sum)/tapply(rowSums(internalData[windowsToUse,3:5]),windows[windowsToUse],sum)
  }),recursive=TRUE)
  avg_posFreq<-cbind(avg_pos,avg_freq)
  avg_posFreq<-t(sapply(sort(avg_posFreq[,1],index.return=T)$ix,function(x) avg_posFreq[x,]))
  avg_minIndex<-which(min(abs(avg_freq-foreground_frequency))==abs(avg_freq-foreground_frequency))

  print(paste("Finding initial peak(s).. min distance in a window of size ",winSize," bp to ",foreground_frequency,": ",min(abs(avg_freq-foreground_frequency)),sep=""))

  for(index in avg_minIndex){
   print(paste("   At (avg(pos) in window): ",round(avg_pos[index])," bp",sep=""))
  }
  
  #order confidence levels
  level<-sort(level)

  res<- identify_peaks(1,length(internalData[,2]),foreground_frequency,level,minWindow,avg_posFreq,bestsize,recurse,forceInclude, allowAdjustment)
  res<-matrix(res[res[,3]<0,],ncol=4)
  ci<-matrix(c(0,0,920,1),nrow=4)
  if(!is.null(dim(res))&&dim(res)[1]>0){
   ci<-matrix(apply(res,1,function(x) t(c(start=ifelse(x[3]<0,internalData[x[1],2],0), stop=ifelse(x[3]<0,internalData[x[1]+x[2]-1,2],0),p.value=ifelse(x[3]<0,-1*(x[3]+x[2]),x[3]),level=x[4] ))),nrow=4)
   print("Found interval:")
#   print(ci)
   for(i in 1:length(ci[1,])){
    print(paste(ci[1,i],"-",ci[2,i]))
   }
  }
#  apply(ci,2,function(x) print(paste(x[1],"-",x[2])))
  list(confidenceInterval=ci,excluded=filtered)
 }else{
  #too few markers
  list(confidenceInterval=t(t(c(0,0,919))),excluded=filtered)
 }
}

filterSamplingv3 <-function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE){
 fs_freqs<-internalData[,3]/rowSums(internalData[,3:5])
 fs_allPos<-internalData[,2]
 fs_n<-length(fs_freqs)
 fs_tested<- rep(FALSE,fs_n)
 fs_filter<- rep(FALSE,fs_n)
 fs_allIndices<-1:fs_n
 fs_chrStart<-min(fs_allPos)
 fs_chrEnd<-max(fs_allPos)

 while(sum(!fs_tested)>0){
  fs_diff<-abs(diff(fs_freqs[!fs_filter]))
  fs_diff<-c(0,fs_diff)+c(fs_diff,0)
  
  fs_curIndex<-fs_allIndices[!fs_filter][!fs_tested[!fs_filter]][which.max(fs_diff[!fs_tested[!fs_filter]])]
  fs_curPos<-internalData[fs_curIndex,2]

#   plot(internalData[!fs_filter,2],fs_diff)
#   abline(v=fs_curPos,col="red")

  #calculate window to use
  fs_start<-max(fs_chrStart,fs_curPos-fs_windowsize/2)
  fs_end<- fs_start + fs_windowsize
  if(fs_end>fs_chrEnd){
   fs_start<-max(fs_chrStart,fs_end-fs_windowsize)
  }
  fs_toUse<-fs_allPos>=fs_start & fs_allPos<=fs_end & fs_allPos != fs_curPos & !fs_filter  
  fs_size<- sum(fs_toUse)
  fs_p<-1
  if(fs_size>3){
   assign("dataset_shoremapmle",internalData[fs_toUse,],".GlobalEnv")
   fs_p.win<- samplefreqs(1,fs_size)
#   fs_p.win<-colSums(fs_data[,3:5])/sum(fs_data[,3:5])
   if(fs_exact){
    sink("/dev/null");
    fs_p<-multinomial.test(c(internalData[fs_curIndex,3:5],recursive=TRUE),prob=fs_p.win)$p.value
    sink();
   }else{
    fs_p1<-fs_p.win[3]
    fs_p2<-fs_p.win[1]/sum(fs_p.win[1:2])
    x<-c(internalData[fs_curIndex,3:5],recursive=TRUE)
    fs_p<-pbinom(x[3]+ifelse(x[3]<sum(x)*fs_p1,1,-1),size=sum(x),prob=fs_p1,lower.tail=x[3]<sum(x)*fs_p1)*pbinom(x[1]+ifelse(x[1]<sum(x[1:2])*fs_p2,1,-1),size=sum(x[1:2]),prob=fs_p2,lower.tail=x[1]<sum(x[1:2])*fs_p2)
   }
  }
  fs_filter[fs_curIndex]<- p.adjust(fs_p,method="holm",n=fs_n)<=fs_limit
  fs_tested[fs_curIndex]<-TRUE
 }
 !fs_filter
}

identify_peaks <- function(indexL,indexH,frequency,level,minWindow,avg_posFreq,bestsize,recurse,forceInclude=TRUE,allowAdjustment=0.05){
 assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
 if(indexH-indexL>max(minWindow,bestsize)){
#  print(paste(indexL,indexH,bestsize,sep=" ### "))
  cur_indices<-indexL:indexH
  avg_pf<-avg_posFreq[avg_posFreq[,1]>=min(dataset_shoremapmle[cur_indices,2]) & avg_posFreq[,1]<=max(dataset_shoremapmle[cur_indices,2]),]
  if(TRUE){ #condition to recect if avg_pf is too distant
   #try to find peaks
   starts<-avg_pf[which(min(abs(avg_pf[,2]-frequency))==abs(avg_pf[,2]-frequency)),1]
   while(length(starts)>0){
    start<- starts[ceiling(length(starts)/2)]
    beststarts<-which(min(abs(dataset_shoremapmle[,2]-start))==abs(dataset_shoremapmle[,2]-start))
    beststart<-beststarts[ceiling(length(beststarts)/2)]
    #see if the frequency needs adjustment
    newFreq<-multll(beststart-floor(bestsize/2),bestsize)$coef[1]
    if(abs(frequency-newFreq)<allowAdjustment){
     print(paste("The goal frequency was adjusted from ",frequency," to ",newFreq,", which is the estimated frequency in the preliminary interval",sep=""))
     frequency<-newFreq
    }
    res<-extend(beststart,bestsize,level[1],frequency,indexL,indexH,minWindow, forceInclude)
    if(length(level)>1){
     for(i in 2:length(level)){
      res2<-c(1,1,618,level[i])
      assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
      if(res[i-1,3]>0){
       res2<-extend(beststart,bestsize,level[i],frequency,indexL,indexH,minWindow, forceInclude)
      }else{
       res2<-extend(res[i-1,1],res[i-1,2],level[i],frequency,indexL,indexH,minWindow, forceInclude)
      }
      res<-rbind(res,res2)
     }
    }
    if(min(res[,3])>0){
     if(length(starts)==1){
      return(matrix(c(1,1,616,1),ncol=4))
      break
     }else{
      starts<-starts[1:length(starts)!=ceiling(length(starts)/2)]
     }
    }else if(recurse){
     if(min(res[,3])<0){
      resL<- identify_peaks(indexL, min(res[res[,3]<0,1]), frequency, level[1], minWindow, avg_posFreq, bestsize, recurse, forceInclude,0.0)
      if(length(level)>1){
       for(i in 2:length(level)){
        res2<-c(1,1,618,level[i])
        assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
        if(res[i-1,3]>0){
         res2<-extend(beststart,bestsize,level[i],frequency,indexL,min(res[res[,3]<0,1]),minWindow, forceInclude)
        }else{
         res2<-extend(resL[i-1,1],resL[i-1,2],level[i],frequency,indexL,min(res[res[,3]<0,1]),minWindow, forceInclude)
        }
        resL<-rbind(resL,res2)
       }
      }
      resH<- identify_peaks(max(res[res[,3]<0,1]+res[res[,3]<0,2]), indexH, frequency, level[1], minWindow, avg_posFreq, bestsize, recurse, forceInclude,0.0)
      if(length(level)>1){
       for(i in 2:length(level)){
        res2<-c(1,1,618,level[i])
        assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
        if(res[i-1,3]>0){
         res2<-extend(beststart,bestsize,level[i],frequency,max(res[res[,3]<0,1]+res[res[,3]<0,2]),indexH,minWindow, forceInclude)
        }else{
         res2<-extend(resH[i-1,1],resH[i-1,2],level[i],frequency,max(res[res[,3]<0,1]+res[res[,3]<0,2]),indexH,minWindow, forceInclude)
        }
        resH<-rbind(resH,res2)
       }
      }

      if(min(resL[,3])<0){
       if(min(resH[,3])<0){
        #both good
        return(rbind(resL,res,resH))
        break
       }else{
        #low good
        return(rbind(resL,res))
        break
       }
      }else if(min(resH[,3])<0){
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
   return(t(as.matrix(c(1,1,616,1))))
  }
 }else{
  #Too small window
  return(t(as.matrix(c(1,1,617,1))))
 }
}

loglikelihood_mult <- function(llm_P1=0.5,llm_err=0.01,llm_index=0,llm_size=0){
 #the loglikelihood function. Returns 110000 for unvalid p-values
# if(P1<0 || P1>1 || err<1e-2 || err>1) {
 if(is.na(llm_P1) || is.na(llm_err) || llm_size<0){
  220000
 } else if(llm_P1<0 || llm_P1>1 || llm_err<0 || llm_err>1) {
  110000
 }else{
  llm_P1<-as.numeric(llm_P1)
  llm_err<-as.numeric(llm_err)
  #print(paste(llm_P1,llm_err,sep="##"))
  llm_p1<- llm_P1*(1-4*llm_err/3)+llm_err/3
  llm_pe<- 2*llm_err/3
  llm_p2<- 1-llm_p1-llm_pe
  llm_p.all <- c(llm_p1,llm_p2,llm_pe)
  -sum(apply(dataset_shoremapmle[llm_index:(llm_index+llm_size-1),],1,function(llm_x){dmultinom(x=c(llm_x[3],llm_x[4],llm_x[5]),prob=llm_p.all,log=TRUE)}))
 }
}

samplefreqs <- function(sf_startPos,sf_size) {
 curIndices<- sf_startPos:(sf_startPos+sf_size-1)
 colSums(dataset_shoremapmle[curIndices,3:5])/sum(dataset_shoremapmle[curIndices,3:5])
}

multll<- function(ml_x,ml_size=10) {
 p.win<-samplefreqs(ml_x,ml_size)
 ml_errEst<-3*p.win[3]/2
 ml_P1est<-(p.win[1]-ml_errEst/3)/(1-4*ml_errEst/3)
 #preventing out of bounds results for non-model cases like a homozygous marker with errors
 ml_errEst<-min(1,max(0,ml_errEst))
 ml_P1est<-min(1,max(0,ml_P1est))
 #calculate likelihood
 ml_min<-loglikelihood_mult(llm_P1=ml_P1est,llm_err=ml_errEst,llm_index=ml_x,llm_size=ml_size)
 #return mle2-like results
 list(coef=c(ml_P1est,ml_errEst),min=ml_min)
}

restrictedModel <- function(P1,x,size) {
 rM_errEst<-0.0001
 if(x>0&&x+size<length(dataset_shoremapmle)){
  rM_curIndices<-x:(x+size-1)
  rM_errEst<-3*sum(dataset_shoremapmle[rM_curIndices,5])/sum(dataset_shoremapmle[rM_curIndices,3:5])/2
 }
 mle2(loglikelihood_mult,optimizer="optimize",start=list(llm_err=rM_errEst),fixed=list(llm_P1=P1,llm_index=x,llm_size=size),lower=0, upper=1)
}

maxConf<-function(x,level=0.95,freq=0,indexL=0,indexH=Inf,minWindow=10,include=c(-1,-1)){
 #function to minimize for the optimization of the interval
 start<-floor(x[1])
 size<-floor(x[2])
# print(paste(start,size,sep=" ## "))
 indexL<- max(1,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 if(size<minWindow){
  140000-size+minWindow
 }else if(start<indexL){
  130000+indexL-start
 }else if(start+size-1>indexH){
  120000-indexH-1+start+size
 }else if(sum(include)>0 && (start>include[1])){
  #given region not included
  110000+start-include[1]
 }else if(sum(include)>0 && (start+size < sum(include))){
  #given region not included
  110000-start-size+sum(include)
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
   if(fit$min>100000){
    res<-fit$min
   }else{
    restrictedFit<- restrictedModel(freq,start,size)
#    p<- pchisq(-2*(fit@min-restrictedFit@min),1)
    p<- pchisq(-2*(fit$min-restrictedFit@min),1)
    if(p<=level){
     res<- -size-p
    }else{
     res<- size+p
    }
   }
   assign("storage_shoremapmle",rbind(storage_shoremapmle,c(start,size,res)),".GlobalEnv")
#   print(paste(start,size,res,sep=" ## "))
   res
  }
 }
}

pInterval<-function(start,size,freq){
 fit<-multll(start,size)
 restrictedFit<- restrictedModel(freq,start,size)
 pchisq(-2*(fit$min-restrictedFit@min),1)
}

extend <- function(beststart,bestsize=10,level=0.99,freq=1,indexL=0,indexH=Inf,minWindow=10,forceInclude=TRUE){
 #given a window it extends this as far as possible to the left and right without exceeding the confidence level
 bestvalue<-Inf
 indexL<- max(1,indexL)
 indexH<- min(length(dataset_shoremapmle[,2]),indexH)
 inclusion<-c(-1,-1)
 if(forceInclude){
  inclusion<-c(beststart,bestsize)
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
 matrix(as.numeric(c(beststart,bestsize,bestvalue,level)),ncol=4)
}

#not better
filterSamplingv2 <- function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE){
 fs_allPos<-internalData$V2
 fs_allIndices<-1:length(fs_allPos)
 fs_chrStart<-min(fs_allPos)
 fs_chrEnd<-max(fs_allPos)
 fs_ret<-sapply(fs_allIndices, function(fs_curIndex){
  fs_curPos<-internalData$V2[fs_curIndex]
  fs_start<-max(fs_chrStart,fs_curPos-fs_windowsize/2)
  fs_end<- fs_start + fs_windowsize
  if(fs_end>fs_chrEnd){
   fs_start<-max(fs_chrStart,fs_end-fs_windowsize)
  }
  fs_data<- internalData[fs_allPos>=fs_start & fs_allPos<=fs_end & fs_allPos != fs_curPos,]
  fs_size<- length(fs_data$V1)
  if(fs_size>3){
   assign("dataset_shoremapmle",fs_data,".GlobalEnv")
   fs_p.win<- samplefreqs(1,fs_size)
   if(fs_exact){
    sink("/dev/null");
    multinomial.test(c(internalData[fs_curIndex,3:5],recursive=TRUE),prob=fs_p.win)$p.value
    sink();
   }else{
    fs_p1<-fs_p.win[3]
    fs_p2<-fs_p.win[1]/sum(fs_p.win[1:2])
    fs_p2Alt<-fs_p.win[3]/sum(fs_p.win[1:2])
    x<-c(internalData[fs_curIndex,3:5],recursive=TRUE)
    pbinom(x[3],size=sum(x),prob=fs_p1)*ifelse(x[1]<x[2],pbinom(x[1],size=sum(x[1:2]),prob=fs_p2),pbinom(x[2],size=sum(x[1:2]),prob=fs_p2Alt))
   }
  }else{
   1
  }
 })
 p.adjust(fs_ret,method="holm")>=fs_limit
}


#not as good
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
#    fs_p2Alt<-fs_p.win[2]/sum(fs_p.win[1:2])
    fs_p.res<-apply(internalData[fs_curIndices,3:5],1,function(x) pbinom(x[3],size=sum(x),prob=fs_p1)*pbinom(x[1],size=sum(x[1:2]),prob=fs_p2,lower.tail=x[1]<sum(x[1:2])*fs_p2))
#    fs_p.res<-apply(internalData[fs_curIndices,3:5],1,function(x) pbinom(x[3],size=sum(x),prob=fs_p1)*ifelse(x[1]<x[2],pbinom(x[1],size=sum(x[1:2]),prob=fs_p2),pbinom(x[2],size=sum(x[1:2]),prob=fs_p2Alt)));
   }
   fs_ret<-c(fs_ret,fs_p.res)
  }else{
   fs_ret<-c(fs_ret,rep(1,fs_size))
  }
 }
 p.adjust(fs_ret,method="holm")>=fs_limit
}

