#returns the reached p-value for the interval along with the positions (start, stop, p.value)
#if the third value is larger than the confidence level, it is not a valid region.
#616 not a valid peak
#617 the examined window was too small compared to the internal windowsize (minWindow)
#919 too few markers after removal of 0-sum-markers, filtering. Too few means that there are less than ten or that next too all are constant

#chromosome,positions, background_count,forground_count and error_count are vectors of the same length
require(bbmle)
#require(EMT)

ShoreMap.confint <- function(chromosome,positions, background_count, foreground_count, error_count, foreground_frequency=1, level=0.99, recurse=FALSE, forceInclude=TRUE, allowAdjustment=0.0, filterOutliers=200000, filterPValue=0.05, winSize=50000, winStep=10000, minMarker=10, minCoverage=0, maxCoverage=Inf, peakFinding=3, peakWinSize=50000, peakWinStep=10000) {
 
 verbose=FALSE

 if(verbose){
  print(paste("foreground_frequency=", foreground_frequency), sep="")
  print(paste("level=", level), sep="")
  print(paste("recurse=", recurse), sep="")
  print(paste("forceInclude=", forceInclude), sep="")
  print(paste("allowAdjustment=", allowAdjustment), sep="")
  print(paste("filterOutliers=", filterOutliers), sep="")
  print(paste("filterPValue=", filterPValue), sep="")
  print(paste("winSize=", winSize), sep="")
  print(paste("winStep=", winStep), sep="")
  print(paste("minMarker=", minMarker), sep="")
  print(paste("minCoverage=", minCoverage), sep="")
  print(paste("maxCoverage=", maxCoverage), sep="")
  print(paste("peakFinding=", peakFinding), sep="")
  print(paste("peakWinSize=", peakWinSize), sep="")
  print(paste("peakWinStep=", peakWinStep), sep="")
 }

# allowAdjustment=0.0
# minMarker=10
# level<-c(0.95,0.99,0.999)
# print(sapply(ls(all.names=TRUE),function(x) eval(parse(text=paste("length(",x,")",sep="")))))
# peakFinding<-3 #3 is boost, 4 is R
# recurse=TRUE
 minMarker<-max(1,minMarker)
 foreground_frequency<-as.numeric(foreground_frequency)
 internalData<- cbind(chromosome,positions,foreground_count,background_count,error_count)
# internalData<- data.frame(V1=chromosome, V2=positions, V3=foreground_count, V4=background_count, V5=error_count)
 internalData<- internalData[rowSums(internalData[,3:5])>minCoverage&rowSums(internalData[,3:5])<maxCoverage,]
 print(paste("Analysing chr ",chromosome[1],", with ",length(chromosome)," (",length(internalData[,1]),") markers for equality to ",foreground_frequency,"(",typeof(foreground_frequency),")",sep=""))
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
 freqs<-internalData[,3]/rowSums(internalData[,3:5])
 assign("i_shoremapmle",0,".GlobalEnv")
 minWindow<-max(minMarker,2)
# bestsize<- ceiling((max(table(sapply(2:length(freqs),function(x) if(freqs[x]==freqs[x-1]){i_shoremapmle}else{assign("i_shoremapmle",i_shoremapmle+1,".GlobalEnv");i_shoremapmle})))+1)/5)*5
 bestsize<-max(minMarker,2)
 bestsize<-max(bestsize,minWindow)
# print(paste("Bestsize:",bestsize))
 if(bestsize<length(dataset_shoremapmle[,1])){
  avg_pos<-c(sapply(seq(0,winSize-1,winStep),function(shift){
   windows<- floor((internalData[,2]+shift)/winSize)
   windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
   tapply(internalData[windowsToUse,2],windows[windowsToUse],mean)
  }),recursive=TRUE)
  avg_freq<-c(sapply(seq(0,winSize-1,winStep),function(shift){
   windows<- floor((internalData[,2]+shift)/winSize)
   windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
   tapply(internalData[windowsToUse,3],windows[windowsToUse],sum)/tapply(rowSums(internalData[windowsToUse,3:5]),windows[windowsToUse],sum)
  }),recursive=TRUE)  
  avg_R<-c(sapply(seq(0,winSize-1,winStep),function(shift){
   windows<- floor((internalData[,2]+shift)/winSize)
   windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
   allele<-tapply(internalData[windowsToUse,3],windows[windowsToUse],sum)
   ref<-tapply(internalData[windowsToUse,4],windows[windowsToUse],sum)
   ret<-pmax(allele/ref,ref/allele)
   ret[ret==1]<-0
   rMax=max(ret[!is.infinite(ret)])
   ret[is.infinite(ret)]<-rMax
   ret
  }),recursive=TRUE)
  avg_boost<-abs(1/(1-max(foreground_frequency,1-foreground_frequency)/pmax(avg_freq,1-avg_freq)))
  boostMax<-max(avg_boost[!is.infinite(avg_boost)])
  avg_boost[is.infinite(avg_boost)]<-boostMax
#  avg_boost<-avg_boost/max(avg_boost)
  avg_posFreq<-cbind(avg_pos,avg_freq,avg_boost,avg_R)
  avg_posFreq<-t(sapply(sort(avg_posFreq[,1],index.return=T)$ix,function(x) avg_posFreq[x,]))
  #avg_minIndex<-which(min(abs(avg_posFreq[,2]-foreground_frequency))==abs(avg_posFreq[,2]-foreground_frequency))





  ci<-matrix(c(0,0,920,1,0),nrow=5)
  if(level[1]<=1){
   
   peak_pos<-c(sapply(seq(0,peakWinSize-1,peakWinStep),function(shift){
    windows<- floor((internalData[,2]+shift)/peakWinSize)
    windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
    tapply(internalData[windowsToUse,2],windows[windowsToUse],mean)
   }),recursive=TRUE)
   peak_freq<-c(sapply(seq(0,peakWinSize-1,peakWinStep),function(shift){
    windows<- floor((internalData[,2]+shift)/peakWinSize)
    windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
    tapply(internalData[windowsToUse,3],windows[windowsToUse],sum)/tapply(rowSums(internalData[windowsToUse,3:5]),windows[windowsToUse],sum)
   }),recursive=TRUE)  
   peak_R<-c(sapply(seq(0,peakWinSize-1,peakWinStep),function(shift){
    windows<- floor((internalData[,2]+shift)/peakWinSize)
    windowsToUse<- windows %in% unique(windows)[table(windows)>minMarker]
    allele<-tapply(internalData[windowsToUse,3],windows[windowsToUse],sum)
    ref<-tapply(internalData[windowsToUse,4],windows[windowsToUse],sum)
    ret<-pmax(allele/ref,ref/allele)
    ret[ret==1]<-0
    rMax=max(ret[!is.infinite(ret)])
    ret[is.infinite(ret)]<-rMax
    ret
   }),recursive=TRUE)
   peak_boost<-abs(1/(1-max(foreground_frequency,1-foreground_frequency)/pmax(peak_freq,1-peak_freq)))
   boostMax<-max(peak_boost[!is.infinite(peak_boost)])
   peak_boost[is.infinite(peak_boost)]<-boostMax

   peak_posFreq<-cbind(peak_pos,peak_freq,peak_boost,peak_R)
   peak_posFreq<-t(sapply(sort(peak_posFreq[,1],index.return=T)$ix,function(x) peak_posFreq[x,]))

   peak_minIndex<-which(peak_posFreq[,peakFinding]==max(peak_posFreq[,peakFinding]))
   print(paste("Finding initial peak(s).. choosen method in a window of size ",peakWinSize," bp with step size of ", peakWinStep, " bp",sep=""))

   for(index in peak_minIndex){
    print(paste("   At (avg(pos) in window): ",round(peak_posFreq[index,1])," bp",sep=""))
   }
  
   #order confidence levels
   level<-sort(level)

   res<- identify_peaks(1,length(internalData[,2]),foreground_frequency,level,minWindow,peak_posFreq[,c(1,peakFinding)],bestsize,recurse,forceInclude, allowAdjustment)
   res<-matrix(res[res[,3]<0,],ncol=4)
   if(!is.null(dim(res))&&dim(res)[1]>0){
    ci<-matrix(apply(res,1,function(x) t(c(start=ifelse(x[3]<0,internalData[max(x[1]-1,1),2]+1,0), stop=ifelse(x[3]<0,internalData[min(x[1]+x[2],length(internalData[,2])),2]-1,0),p.value=ifelse(x[3]<0,x[4]+x[3]+x[2],x[3]),level=x[4],nrOfMarkers= x[2]))),nrow=5)
    print("Found interval:")
#    print(ci)
    for(i in 1:length(ci[1,])){
     print(paste(ci[1,i],"-",ci[2,i],"level:",ci[4,i]))
    }
   }
  }
#  plot(ci)
#  apply(ci,2,function(x) print(paste(x[1],"-",x[2])))
  list(confidenceInterval=ci,excluded=filtered,averaged=avg_posFreq)
 }else{
  #too few markers
  list(confidenceInterval=matrix(c(0,0,919,1,0),nrow=5),excluded=filtered,averaged=c(-1,-1))
 }
}

#version5
filterSampling<-function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE){
 size<-length(internalData[,2])
 chrStart<-min(internalData[,2])
 chrEnd<-max(internalData[,2])
 indices<-1:size
 largeWindow<-fs_windowsize*4 #use the double size to make sure the window is enclosed
 windows1<-floor((internalData[,2])/largeWindow)
 uniqueWin1<-unique(windows1)
 windows2<-floor((internalData[,2]+largeWindow/2)/largeWindow)
 uniqueWin2<-unique(windows2)
 freqs<-internalData[,3]/rowSums(internalData[,3:5])
 diffDataRaw<-abs(diff(freqs))
 diffDataRaw<-c(diffDataRaw[1],diffDataRaw)+c(diffDataRaw,diffDataRaw[size-1])
 diffDataMod<-diffDataRaw
 allPos<-internalData[,2]

 #corresponds to adjusting the p-values below with respect to size nr of tests
 limit<-fs_limit/size

 #use two frames and choose the frame closest to the current marker position
 data1<-tapply(indices,windows1,function(x) data.frame(internalData[x,2] ,internalData[x,3] ,internalData[x,4] ,internalData[x,5] ,rep(FALSE,length(x))))
 data2<-tapply(indices,windows2,function(x) data.frame(internalData[x,2] ,internalData[x,3] ,internalData[x,4] ,internalData[x,5] ,rep(FALSE,length(x))))
 filtered<-c()
 for(i in indices){
  #get marker to test this iteration
  curIndex<-which.max(diffDataMod)

  curPos<-allPos[curIndex]
  if(curPos==19806951) break
  #start and end of window
  start<-max(chrStart,curPos-fs_windowsize/2) #0
  end<- start + fs_windowsize #0
  if(end>chrEnd){ #0
   start<-max(chrStart,end-fs_windowsize)
  }

  #decide which frame and which windows to use
  curWin1<-curPos/largeWindow
  curWin2<-curWin1+0.5
  use1<-abs((curWin1%%1)-0.5)<abs((curWin2%%1)-0.5)
  curWin1<-which(uniqueWin1==floor(curWin1))
 # if(curWin1==178) break
  curWin2<-which(uniqueWin2==floor(curWin2))
#  print(paste(curWin1,curWin2))
  red<-c()
  curSize<-0
  if(use1){
   #include markers within window
   toUse<-data1[[curWin1]][,1]>=start &data1[[curWin1]][,1]<=end
   #add four closest markers
   toUse<- toUse | 1:length(toUse) %in% sort(abs(data2[[curWin2]][,1]-curPos),index.return=TRUE)$ix[1:min(3,length(toUse))]
   curSize<-sum(toUse)
   red<-data1[[curWin1]][toUse,]
  }else{
   #include markers within window
   toUse<-data2[[curWin2]][,1]>=start &data2[[curWin2]][,1]<=end
   #add four closest markers
   toUse<- toUse | 1:length(toUse) %in% sort(abs(data2[[curWin2]][,1]-curPos),index.return=TRUE)$ix[1:min(3,length(toUse))]
   curSize<-sum(toUse)
   red<-data2[[curWin2]][toUse,]
  }
  p<-if(curSize>3){
   x<-c(red[red[,1]==curPos,2:4],recursive=TRUE)
   red2<- red[red[,1]!=curPos & !red[,5],]
   if(nrow(red2)>0){
    p.win<-colSums(red2[,2:4])
    p.win<-p.win+1/500 
    p.win<-p.win/sum(p.win) #0

    
#   if(fs_exact){
#    sink("/dev/null");
#    p<-multinomial.test(x,prob=fs_p.win)$p.value
#    sink();
#   }else{
    fs_p1<-p.win[3] #0
    fs_p2<-p.win[1]/sum(p.win[1:2]) #0
    pbinom(x[3]+ifelse(x[3]<sum(x)*fs_p1,1,-1),size=sum(x),prob=fs_p1,lower.tail=x[3]<sum(x)*fs_p1)*pbinom(x[1]+ifelse(x[1]<sum(x[1:2])*fs_p2,1,-1),size=sum(x[1:2]),prob=fs_p2,lower.tail=x[1]<sum(x[1:2])*fs_p2) #0.001
#   }
   }else{
    1
   }
  }else{
   1
  }
  #judgement
  if(p<=limit){ #0.011

   #mark outlier
   if (length(data1[[curWin1]][data1[[curWin1]][,1]==curPos,5]) != 0) {
    data1[[curWin1]][data1[[curWin1]][,1]==curPos,5]<-TRUE #0.001
   }
   if (length(data2[[curWin2]][data2[[curWin2]][,1]==curPos,5]) != 0) {
    data2[[curWin2]][data2[[curWin2]][,1]==curPos,5]<-TRUE #0.001
   }
   diffDataMod[curIndex]<--2 #0.001
  
   #recalculate diff values for neighboring markers
   before<-curIndex-1 #0
   while(sum(before==filtered)>0){ #~0 with 57 markers in filtered
    before<-before-1
   }
   after<-curIndex+1 #0
   while(sum(after==filtered)>0){ #0 with no markers in filtered
    after<-after+1
   }

   if(before<1){ #0
    #first marker
    diffDataRaw[after]<-2*(diffDataRaw[after]-abs(freqs[after]-freqs[curIndex]))
    if(diffDataMod[after]>=0){
     diffDataMod[after]<-diffDataRaw[after]
    }
   }else if(after>size){ #0
    #last marker
    diffDataRaw[before]<-2*(diffDataRaw[before]-abs(freqs[before]-freqs[curIndex]))
    if(diffDataMod[before]>=0){
     diffDataMod[before]<-diffDataRaw[before]
    }
   }else{
    flank<-abs(freqs[before]-freqs[after]) #0

    diffDataRaw[before]<-diffDataRaw[before]-abs(freqs[before]-freqs[curIndex])+flank #0
    diffDataRaw[after]<-diffDataRaw[after]-abs(freqs[after]-freqs[curIndex])+flank
    if(diffDataMod[after]>=0){
     diffDataMod[after]<-diffDataRaw[after]
    }
    if(diffDataMod[before]>=0){
     diffDataMod[before]<-diffDataRaw[before]
    }
   }

   #add the position to the filtered positions
   filtered<-c(filtered,curIndex) #0
  }else{
   diffDataMod[curIndex]<--1
  } 
 }
 diffDataMod!=-2
}


identify_peaks <- function(indexL,indexH,frequency,level,minWindow,avg_posFreq,bestsize,recurse,forceInclude=TRUE,allowAdjustment=0.05){
 assign("storage_shoremapmle",matrix(c(-1,-1,-1),nrow=1),".GlobalEnv")
 if(indexH-indexL>min(minWindow,bestsize)){ #too small window
#  print(paste(indexL,indexH,bestsize,sep=" ### "))
  cur_indices<-indexL:indexH
  avg_toUse<-avg_posFreq[,1]>=min(dataset_shoremapmle[cur_indices,2]) & avg_posFreq[,1]<=max(dataset_shoremapmle[cur_indices,2])
  if(sum(avg_toUse)>1){ #too few windowed markers
   avg_pf<-avg_posFreq[avg_toUse,]
   #try to find peaks
#   starts<-avg_pf[which(min(abs(avg_pf[,2]-frequency))==abs(avg_pf[,2]-frequency)),1]
   starts<-avg_pf[which(avg_pf[,2]==max(avg_pf[,2])),1]
   while(length(starts)>0){
    start<- starts[ceiling(length(starts)/2)]
    beststarts<-which(min(abs(dataset_shoremapmle[,2]-start))==abs(dataset_shoremapmle[,2]-start))
    beststart<-beststarts[ceiling(length(beststarts)/2)]
    beststart<-min(length(dataset_shoremapmle[,1])-bestsize,beststart)
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
   #Too few windowed markers
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
     res<- -size-(level-p)
    }else{
     res<- size+(level-p)
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
  
  while(i<minWindow*2&&curvalue<lastvalue){
   lastvalue<-curvalue
   i<-i+1
   curvalue<-maxConf(c(beststart-floor(i/2),bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
  }
  i<-i-1
  beststart<-beststart-floor(i/2)
  bestsize<-bestsize+i
  if(i%%2==0 && i<2*minWindow){
   #if the extension breaks down on the right, continue on to expand left
   i<-1
   lastvalue<-maxConf(c(beststart,bestsize),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   curvalue<-maxConf(c(beststart-i,bestsize+i),level=level,freq=freq,indexL=indexL,indexH=indexH,minWindow=minWindow,include=inclusion)
   while(i<minWindow&&curvalue<lastvalue){
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
   while(i<minWindow&&curvalue<lastvalue){
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




