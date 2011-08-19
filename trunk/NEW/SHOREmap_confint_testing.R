winSize<-5e5
winStep<-1e4
minMarkersInWindow<-20
sapply(seq(50000,1000000,50000),function(winSize){
shifts<-seq(0,winSize-1,winStep)

unsorted<-apply(cbind(shifts,shifts),1,function(x){
shift<-x[1]
#filter??
passedData<-internalData
windows<-floor((internalData[,2]+shift)/winSize)
#remove small windows
passedCount<- table(windows)
windowsToUse<-as.numeric(rownames(passedCount)[passedCount>minMarkersInWindow])
markersToUse<- windows %in% windowsToUse
passedData<-passedData[markersToUse,]
windows<-windows[markersToUse]
boost<-abs(1/(1-foreground_frequency/sapply(tapply(passedData[,3],windows,sum)/tapply(rowSums(passedData[,3:5]),windows,sum),function(f) max(f,1-f))))
pos<-tapply(passedData[,2],windows,mean)
#matrix(pos,boost,ncol=2)
cbind(pos,boost)
})
m<-t(sapply(unsorted,function(x) x[which.max(x[,2]),]))
cbind(winSize,m[m[,2]==max(m[,2]),1])
})
freqs<-internalData[,3]/rowSums(internalData[,3:5])
plot(internalData[,2],freqs,xlim=c(2.5e6,3.5e6))

pos<-c(sapply(unsorted,function(i) i[,1]),recursive=T)
boost<-c(sapply(unsorted,function(i) i[,2]),recursive=T)
# aa<-c(sapply(unsorted,function(i) i[,3]),recursive=T) #needed for calculation of allele
o<-sort.int(pos,index.return=T)
# sorted<-t(sapply(o$ix,function(i) c(pos[i],tp[i],aa[i]))) #needed for calculation of allele
sorted<-t(sapply(o$ix,function(i) c(pos[i],boost[i])))
sorted[which.max(sorted[,2]),]


plot(internalData[,2],freqs)
plot(internalData[,2],freqs,xlim=c(2.5e6,3.5e6))
lines(sorted[,1],sorted[,2]/max(sorted[,2]),col="red")


pscore.loess<-loess(pscore~x,span=0.1,data.frame(x=sorted[,1],pscore=sorted[,2]))
pscore.predict<- predict(pscore.loess,data.frame(x=sorted[,1]))
sorted[which.max(pscore.predict),]



lines(sorted[,1],pscore.predict/max(pscore.predict),col="orange",type="l")
