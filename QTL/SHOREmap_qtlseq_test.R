



#linear fitting

peakFun<-function(x,k1,m,b,k2){
 ifelse(x<b,k1*x+m,k2*x+(k1-k2)*b+m)
}

platauFun<-function(x,k1,m,b1,b2,k2){
 ifelse(x<b1,k1*x+m,ifelse(x<b2,k1*b1+m,k2*(x-b2)+k1*b1+m))
}




for (chr in 1:5){
png(paste("summary_chr",chr,".png",sep=""))
x<-which(data[,1]==chr)
y.loess<-loess(y~x,span=0.75,data.frame(x=data[x,2],y=p[x]))
y.predict<-predict(y.loess,data.frame(x=data[x,2]))


#plot(data[x,2],data[x,3],col="grey",type="l",yaxt="n",ylab="")
#axis(4)
#mtext("frequency",side=4,line=3)
#lines(data[x,2],data[x,7],col="black",type="l")
#par(new=T)
#plot(data[x,2],p[x],xlab="pos",ylab="log(p)",main=paste("chr",chr))
plot(data[x,2],y.predict,col="red",type="l",xlab="pos",ylab="p smoothed",main=paste("chr",chr))
#lines(data[x,2],y.predict,col="red",type="l")
for(pos in qtls[qtls$chr==chr,1]){
 abline(v=pos,col="green")
}
dev.off()

}

sapply(1:length(data[,1]),mult_one_fun,r)

assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1,-1),nrow=1),".GlobalEnv")
#assign("shoremap_qtsData",data,".GlobalEnv")
#  r<- 250000*100/nrOfPlants

reduceDataForChr<-function(data,winSize,chr){
 windows<-floor(data[,2]/winSize)
 redData<-apply(data[data[,1]==chr,c(4:6,8:10)],2,function(x) tapply(x,windows,sum))
 cbind(rep(chr,length(unique(windows))),unique(windows)*winSize+winSize/2,apply(redData,1,function(x) x[1]/sum(x[1:3])),redData[,1:3],apply(redData,1,function(x) x[4]/sum(x[4:6])),redData[,4:6])
}


mult_two_fun <- function(locus,r,limit=0.9,memmory=3){
 if((locus%%1000)==0){
  print(locus)
 }
 dist<-abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2])
 assign("shoremap_q",sapply(dist,q_fun,r),".GlobalEnv")
# assign("shoremap_mod",ppois(dist,r,lower.tail=F),".GlobalEnv")
 assign("shoremap_toUse",shoremap_q>limit,".GlobalEnv")
# assign("shoremap_toUse",ppois(dist,r,lower.tail=F)>0,".GlobalEnv")
 if(sum(shoremap_toUse)>1){
  mle2(ll_two_fun,optimizer="optimize",start=list(f=0.5),fixed=list(locus=locus,r=r,memmory=memmory),lower=0,upper=1)@coef
 }else{
  sum(shoremap_qtlData[locus,c(4,8)])/sum(shoremap_qtlData[locus,c(4,5,8,9)])
 }

}

mult_one_fun<-function(locus,r,limit=0.9,memmory=1){
 if((locus%%1000)==0){
  print(locus)
 }
 cols=c(4,5)
 if(memmory==2){
  cols=c(8,9)
 }
 dist<-abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2])
 assign("shoremap_q",sapply(dist,q_fun,r),".GlobalEnv")
# assign("shoremap_mod",ppois(dist,r,lower.tail=F),".GlobalEnv")
 assign("shoremap_toUse",shoremap_q>limit,".GlobalEnv")
# assign("shoremap_toUse",ppois(dist,r,lower.tail=F)>0,".GlobalEnv")
 if(sum(shoremap_toUse)>1){
  mle2(ll_one_fun,optimizer="optimize",start=list(f=0.5),fixed=list(locus=locus,r=r,memmory=memmory),lower=0,upper=1)@coef
 }else{
  shoremap_qtlData[locus,cols[1]]/sum(shoremap_qtlData[locus,cols])
 }
}

maxConf<-function(x,minDiff,level,validMin=1,validMax=Inf,maxStart=Inf,minStop=-Inf){
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
 }else{
  interval<-start:stop
  f1<-estimateFreq_one(interval,1)
  f2<-estimateFreq_one(interval,2)
  diff<-abs(f1-f2)
  if(diff>=minDiff){
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
    print(paste(p,start,stop,shoremap_qtlData[start,2],shoremap_qtlData[stop,2],ll,sep=" ## "))
    ll
   }else{
    shoremap_qtlData[stop,2]-shoremap_qtlData[start,2]
   }
  }
 }
}

ll_rest_two<-function(f1,f2,minDiff){
 if(f1<0||f1>1){
  110000+abs(f1-0.5)
 }else if(f2<0 || f2>1){
  120000+abs(f2-0.5)
 }else if(abs(f1-f2)<minDiff){
  130000+abs(f1-f2)
 }else{
  ll_one_fun(f1,shoremap_interval,1)+ll_one_fun(f2,shoremap_interval,2)
 }
}

estimateFreq_one<-function(interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 f<-colSums(shoremap_qtlData[interval,cols])
 f[1]/sum(f)
}

estimateFreq_two<-function(interval){
 cols<-c(4,5,8,9)
 f<-colSums(shoremap_qtlData[interval,cols])
 (f[1]+f[3])/sum(f)
}


ll_one_fun<-function(f,interval,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
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

ll_two_fun<-function(f,interval,memory=3){
 cols1<-c(4,5)
 cols2<-c(8,9)
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else{
  l<- ll_one_fun(fc,interval,1)+ll_one_fun(fc,interval,2)
  l
 }
}

ll_two_fun_withQ<-function(f,locus,r=50000,memory=3){
 cols1<-c(4,5)
 cols2<-c(8,9)
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory,4]
 }else{
  l<- ll_one_fun_withQ(fc,locus,r,1)+ll_one_fun_withQ(fc,locus,r,2)
#  print(paste(fc,l,sep=" ## "))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,locus,memmory,l)),".GlobalEnv")
  l
 }
}

ll_one_fun_withQ<-function(f,locus,r=50000,memory=1){
 cols<-c(4,5)
 if(memory==2){
  cols<-c(8,9)
 }
# print(summary(q))
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memory,4]
 }else{
#  r<- 250000*100/nrOfPlants
#  q<- sapply(abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2]),q_fun,r)
  l<- -sum(apply(cbind(shoremap_q[shoremap_toUse],shoremap_qtlData[shoremap_toUse,cols]),1,function(x) pa_fun(fc,x[1],x[2],x[3])))
#  print(paste(fc,l,sep=" ## "))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,locus,memmory,l)),".GlobalEnv")
  l
 }
}

pa_fun<-function(f,q,p1,p2){
 #mod: supposed to be the probability that less than
 res<-0 
 if(q==0){
  #0
  res<-(1-f)^p1*f^p2
 }else{
  res<-(q*f^p1*(1-f)^p2+(1-q)*(1-f)^p1*f^p2)
 }
 ifelse(res>0,log(res),-61600)
}

q_fun<-function(d,r){
 exp(-d/r)
}

















plot(data[,2],data[,3])
lines(data[,2],data[,7],col="red",type="l")
lines(pos,freqs[,1],col="green",type="l")
lines(pos,freqs[,2],col="blue",type="l")

