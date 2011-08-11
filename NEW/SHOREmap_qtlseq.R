args<-commandArgs(trailingOnly=TRUE)

high_file<-args[1]
low_file<-args[2]
qtl_file<-args[3]

#prep data
hs<-read.table("high.seq")
ls<-read.table("low.seq")

##merge
hs.freq<- apply(hs,1,function(x) x[3]/sum(x[3:5]))
ls.freq<- apply(ls,1,function(x) x[3]/sum(x[3:5]))
data<- cbind(hs[,1:2],hs.freq,hs[,3:5],ls.freq,ls[,3:5])
#fisher test
p<-apply(data,1,function(x) fisher.test(matrix(as.numeric(x[c(4,5,8,9)]),nrow=2))$p.value)
lp<- -log(p)

qtls<-data.frame(pos=c(23104253,13746533,15359142,10500440,25490969),chr=c(1,2,4,5,5))
pdata<-cbind(data[,1:2],p,data[,3]-data[,7])

winSize<-500000
winStep<-50000
minMarkersInWindow<-20

for (chr in 1:5){
 x<-which(data[,1]==chr)
 shifts<-seq(0,winSize-1,winStep)

 unsorted<-sapply(shifts,function(shift){
  markerCount<-table(floor((pdata[x,2]+shift)/winSize))
  markerCount<-cbind(as.numeric(rownames(markerCount)),markerCount)

  curData<-pdata[x,][pdata[x,3]<0.2,]
  windows<-floor((curData[,2]+shift) /winSize)
  passedCount<- table(windows)
  windowsToUse<-as.numeric(rownames(passedCount)[passedCount>minMarkersInWindow])
  markersToUse<- windows %in% windowsToUse

  curData<-curData[markersToUse,]
  windows<-windows[markersToUse]

  aa<-tapply(curData[,4],windows,function(window){sum(window)})/markerCount[markerCount[,1] %in% unique(windows),2]
#  tp<-tapply(curData[,3],windows,function(window){(sum(window<0.2)^2)/(length(pdata[x,3])*sum(window[window<0.2]))})
  tp<-tapply(curData[,3],windows,function(window){(sum(window<0.2)^2)/(sum(window[window<0.2]))})/markerCount[markerCount[,1] %in% unique(windows),2]
#  pos<-(unique(windows)*winSize)+(winSize+shift)/2
  pos<-tapply(curData[,2],windows,mean)
  cbind(pos,tp,aa)
 },simplify=F)
 pos<-c(sapply(unsorted,function(i) i[,1]),recursive=T)
 tp<-c(sapply(unsorted,function(i) i[,2]),recursive=T)
 aa<-c(sapply(unsorted,function(i) i[,3]),recursive=T)
 o<-sort.int(pos,index.return=T)
 sorted<-t(sapply(o$ix,function(i) c(pos[i],tp[i],aa[i])))

# y.loess<-loess(y~x,span=0.75,data.frame(x=sorted[,1],y=sorted[,2]))
# y.predict<- predict(y.loess,data.frame(x=sorted[,1]))


 print(paste("chr:",chr,"no Smooth"))
 print(sorted[which.max(sorted[,2]),1])
# print(paste(sorted[which(diff(sign(diff(y.predict)))==-2)+1,1]))

 png(paste("pointypeak_chr",chr,"_smooth_avgPos_winsize",winSize,"_winstep",winStep,".png",sep=""))
# plot(sorted[,1],sorted[,2]*sorted[,3],type="l",col="green",xlab="",ylab="",yaxt="n",xaxt="n")
# mtext("allele score",side=4,line=3)
# axis(4)
# par(new=T)
 plot(sorted[,1],sorted[,2],type="l",main=paste("pointy peak, chr",chr,"(winsize: ",winSize," bp, winstep: ",winStep," bp)"),xlab="pos",ylab="allele score")

# plot(sorted[,1],sorted[,2],type="l",main=paste("pointy peak, chr",chr,"(winsize: ",winSize," bp, winstep: ",winStep," bp)"),xlab="pos",ylab="p score")
# lines(sorted[,1],y.predict,col="green",type="l")

# y.loess<-loess(y~x,span=0.5,data.frame(x=sorted[,1],y=sorted[,2]))
# y.predict<- predict(y.loess,data.frame(x=sorted[,1]))
# lines(sorted[,1],y.predict,col="red",type="l")
# print(paste("chr:",chr,"span:",0.5))
# print(paste(sorted[which(diff(sign(diff(y.predict)))==-2)+1,1]))

 y.loess<-loess(y~x,span=0.6,data.frame(x=sorted[,1],y=sorted[,2]))
 y.predict<- predict(y.loess,data.frame(x=sorted[,1]))
 lines(sorted[,1],y.predict,col="green",type="l")
 maxIndex<-which(diff(sign(diff(c(-Inf,y.predict,-Inf))))==-2)
 minIndex<-which(diff(sign(diff(c(-Inf,-y.predict,-Inf))))==-2)
 
 minShift<-0
 if(maxIndex[1]<minIndex[1]){
  #starts with a local maxima
  minShift<--1
 }else{
  #starts with a local minima
  minShift<-0
 }
 for(i in 1:length(maxIndex)){
  abline(v=sorted[maxIndex[i],1],col="blue")
  if(y.predict[maxIndex[i]]*100>max(y.predict[maxIndex])){
   checkLow<-FALSE
   checkHigh<-FALSE
   minI<-i+minShift
   if(minI>0&minI<=length(minIndex)){
    checkLow<-y.predict[maxIndex[i]]>2*y.predict[minIndex[minI]]
   }else{
    checkLow<-TRUE
   }
   minI<-minI+1
   if(minI>0&minI<=length(minIndex)){
    checkHigh<-y.predict[maxIndex[i]]>2*y.predict[minIndex[minI]]
   }else{
    checkHigh<-TRUE
   }
   if(checkLow && checkHigh){
    points(sorted[maxIndex[i],1],y.predict[maxIndex[i]],col="blue",pch=8,lwd=3,cex=2)
   }
  }
 }

 print(paste("chr:",chr,"span:",0.6))
 print(paste(sorted[which(diff(sign(diff(y.predict)))==-2)+1,1]))

# y.loess<-loess(y~x,span=0.9,data.frame(x=sorted[,1],y=sorted[,2]))
# y.predict<- predict(y.loess,data.frame(x=sorted[,1]))
# lines(sorted[,1],y.predict,col="blue",type="l")
# print(paste("chr:",chr,"span:",0.9))
# print(paste(sorted[which(diff(sign(diff(y.predict)))==-2)+1,1]))



 for(pos in qtls[qtls$chr==chr,1]){
  abline(v=pos,col="red")
 }
# par(new=T)

 dev.off()

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

assign("shoremap_qtlmem",matrix(c(-1,-1,-1,-1),nrow=1),".GlobalEnv")
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

ll_two_fun<-function(f,locus,r=50000,memmory=3){
 cols1<-c(4,5)
 cols2<-c(8,9)
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memmory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memmory,4]
 }else{
  l<- ll_one_fun(fc,locus,r,1)+ll_one_fun(fc,locus,r,2)
#  print(paste(fc,l,sep=" ## "))
  assign("shoremap_qtlmem",rbind(shoremap_qtlmem,c(fc,locus,memmory,l)),".GlobalEnv")
  l
 }
}

ll_one_fun<-function(f,locus,r=50000,memmory=1){
 cols<-c(4,5)
 if(memmory==2){
  cols<-c(8,9)
 }
# print(summary(q))
 nrOfPlants<-500
 fc<-floor(f*nrOfPlants)/nrOfPlants
 if(fc<0 || fc>1){
  110000+abs(fc-0.5)
 }else if(sum(shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memmory )==1) {
  shoremap_qtlmem[shoremap_qtlmem[,1]==fc &shoremap_qtlmem[,2]==locus &shoremap_qtlmem[,3]==memmory,4]
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

