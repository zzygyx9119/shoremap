#returns the reached p-value for the interval along with the positions (start, stop, p.value)
#if the third value is larger than the confidence level, it is not a valid region.
#616 not a valid peak
#617 the examined window was too small compared to the internal windowsize (minWindow)
#919 too few markers after removal of 0-sum-markers, filtering. Too few means that there are less than ten or that next too all are constant

#chromosome,positions, background_count,forground_count and error_count are vectors of the same length
require(bbmle)
require(EMT)

ShoreMap.confint <- function(chromosome,positions, background_count, foreground_count, error_count, foreground_frequency=1, level=0.99, recurse=FALSE, forceInclude=TRUE, allowAdjustment=0.0, filterOutliers=200000, filterPValue=0.05, winSize=50000, winStep=10000, minMarker=10, minCoverage=0, peakFinding=3, peakWinSize=50000, peakWinStep=10000) {
# allowAdjustment=0.0
# minMarker=10
# level<-c(0.95,0.99,0.999)
# print(sapply(ls(all.names=TRUE),function(x) eval(parse(text=paste("length(",x,")",sep="")))))
# peakFinding<-3 #3 is boost, 4 is R
 
 minMarker<-max(1,minMarker)
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
  avg_boost<-abs(1/(1-foreground_frequency/pmax(avg_freq,1-avg_freq)))
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
   peak_boost<-abs(1/(1-foreground_frequency/pmax(peak_freq,1-peak_freq)))
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
    ci<-matrix(apply(res,1,function(x) t(c(start=ifelse(x[3]<0,internalData[max(x[1]-1,1),2]+1,0), stop=ifelse(x[3]<0,internalData[min(x[1]+x[2],length(internalData[,2])),2]-1,0),p.value=ifelse(x[3]<0,-1*(x[3]+x[2]),x[3]),level=x[4],nrOfMarkers= x[2]))),nrow=5)
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

#version4
filterSamplingLarge<-function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE){
 size<-length(internalData[,2])
 indices<-1:size
 largeWindow<-fs_windowsize*8
 shift<-0
 windows1<-floor((internalData[,2]+shift)/largeWindow)
 shift<-fs_windowsize*4
 windows2<-floor((internalData[,2]+shift)/largeWindow)

 c(tapply(indices,windows1,function(x) filterSampling_sub(internalData[x,],fs_windowsize,fs_limit,fs_exact,size)),recursive=TRUE) & c(tapply(indices,windows2,function(x) filterSampling_sub(internalData[x,],fs_windowsize,fs_limit,fs_exact,size)),recursive=TRUE)
}

filterSampling <-function(internalData,fs_windowsize=200000,fs_limit=0.05,fs_exact=FALSE,fs_totSize){
 fs_freqs<-internalData[,3]/rowSums(internalData[,3:5])
 #fs_allPos<-internalData[,2]
 fs_n<-length(fs_freqs)
 fs_tested<- rep(FALSE,fs_n)
 #fs_filter<- rep(FALSE,fs_n)
 fs_allIndices<-1:fs_n

 fs_chrStart<-min(internalData[,2])
 fs_chrEnd<-max(internalData[,2])
 fs_diff<-abs(diff(fs_freqs))
 fs_diff<-c(fs_diff[1],fs_diff)+c(fs_diff,fs_diff[fs_n-1])
  

 fs_data<-data.frame(fs_diff,fs_freqs,fs_tested,internalData[,2:5],rep(FALSE,fs_n))
 fs_notFiltered<-1:fs_n
 fs_filtered<-c()
 min_diff<-Inf
# system.time({
 while(sum(!fs_data[,3])>0){

#  system.time({
  fs_curIndex<-which.max(fs_data[,1]) #0.001
#  fs_curPos<-fs_data[fs_curIndex,4]
  fs_curPosData<-c(fs_data[fs_curIndex,],recursive=TRUE) #0.002

#   plot(internalData[!fs_filter,2],fs_diff)
#   abline(v=fs_curPos,col="red")

  #calculate window to use
  fs_start<-max(fs_chrStart,fs_curPosData[4]-fs_windowsize/2) #0
  fs_end<- fs_start + fs_windowsize #0
  if(fs_end>fs_chrEnd){ #0
   fs_start<-max(fs_chrStart,fs_end-fs_windowsize)
  }
  fs_toUse<-fs_data[,4]>=fs_start & fs_data[,4]<=fs_end & !fs_data[,8] & fs_data[,4]!=fs_curPos #0.021

  fs_p<-1 #0
  if(sum(fs_toUse)>3){ #0.001
#   fs_red<-fs_data[fs_toUse,] #0.047
#   fs_red<-fs_red[fs_red[,4]!=fs_curPos & !fs_red[,8],] #0.001

   fs_p.win<-colSums(fs_data[fs_toUse,5:7]) #0.022
   fs_p.win<-fs_p.win/sum(fs_p.win) #0
#   fs_p.win<-colSums(fs_data[,3:5])/sum(fs_data[,3:5])
   if(fs_exact){
    sink("/dev/null");
    fs_p<-multinomial.test(c(internalData[fs_curIndex,3:5],recursive=TRUE),prob=fs_p.win)$p.value
    sink();
   }else{
    fs_p1<-fs_p.win[3] #0
    fs_p2<-fs_p.win[1]/sum(fs_p.win[1:2]) #0
    x<-fs_curPosData[5:7] #0
    fs_p<-pbinom(x[3]+ifelse(x[3]<sum(x)*fs_p1,1,-1),size=sum(x),prob=fs_p1,lower.tail=x[3]<sum(x)*fs_p1)*pbinom(x[1]+ifelse(x[1]<sum(x[1:2])*fs_p2,1,-1),size=sum(x[1:2]),prob=fs_p2,lower.tail=x[1]<sum(x[1:2])*fs_p2) #0.001
   }
  }
#  })

  if(p.adjust(fs_p,method="holm",n=fs_totSize)<=fs_limit){ #0.011
#   if(fs_data[fs_curIndex,1]<min_diff){
#    min_diff<-fs_data[fs_curIndex,1]
#    print(fs_data[fs_curIndex,])
#   }
   #mark outlier
   fs_data[fs_curIndex,c(1,3,8)]<-c(-1,TRUE,TRUE) #0.010
  
   #recalculate diff values for neighboring markers
   fs_before<-fs_curIndex-1 #0
   while(sum(fs_before==fs_filtered)>0){ #~0 with no markers in fs_filtered
    fs_before<-fs_before-1
   }
   fs_after<-fs_curIndex+1 #0
   while(sum(fs_after==fs_filtered)>0){ #0 with no markers in fs_filtered
    fs_after<-fs_after+1
   }
   if(fs_before<1){
    #first marker
    fs_data[fs_after,1]<-2*(fs_data[fs_after,1]-abs(fs_data[fs_after,2]-fs_curPosData[2]))
   }else if(fs_after>fs_n){
    #last marker
    fs_data[fs_before,1]<-2*(fs_data[fs_before,1]-abs(fs_data[fs_before,2]-fs_curPosData[2]))
   }else{
    fs_flank<-abs(fs_data[fs_before,2]-fs_data[fs_after,2]) #0.001

    fs_data[fs_before,1]<-fs_data[fs_before,1]-abs(fs_data[fs_before,2]-fs_curPosData[2])+fs_flank #0.011

    fs_data[fs_after,1]<-fs_data[fs_after,1]-abs(fs_data[fs_after,2]-fs_curPosData[2])+fs_flank #0.012
   }



   #add the position to the filtered positions
   fs_filtered<-c(fs_filtered,fs_curIndex) #0
  }else{
   fs_data[fs_curIndex,c(1,3)]<-c(-1,TRUE) #0.011
  }
#  sum(!fs_data[,3])>0
 }
# })
 !fs_data[,8]
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
#  system.time({
  fs_diff<-abs(diff(fs_freqs[!fs_filter])) #0.038
  fs_diff<-c(0,fs_diff)+c(fs_diff,0) #0.043
  
  fs_curIndex<-fs_allIndices[!fs_filter][!fs_tested[!fs_filter]][which.max(fs_diff[!fs_tested[!fs_filter]])] #0.037
  fs_curPos<-internalData[fs_curIndex,2]

#   plot(internalData[!fs_filter,2],fs_diff)
#   abline(v=fs_curPos,col="red")

  #calculate window to use
  fs_start<-max(fs_chrStart,fs_curPos-fs_windowsize/2)
  fs_end<- fs_start + fs_windowsize
  if(fs_end>fs_chrEnd){
   fs_start<-max(fs_chrStart,fs_end-fs_windowsize)
  }
  fs_toUse<-fs_allPos>=fs_start & fs_allPos<=fs_end & fs_allPos != fs_curPos & !fs_filter  #0.021
  fs_size<- sum(fs_toUse)
  fs_p<-1
  if(fs_size>3){
   assign("dataset_shoremapmle",internalData[fs_toUse,],".GlobalEnv") #0.031
   fs_p.win<- samplefreqs(1,fs_size) #0.003
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
  fs_filter[fs_curIndex]<- p.adjust(fs_p,method="holm",n=fs_n)<=fs_limit #0.011
  fs_tested[fs_curIndex]<-TRUE #0
#  })
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


