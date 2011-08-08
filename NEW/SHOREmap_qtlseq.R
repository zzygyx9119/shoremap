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

qtls<-data.frame(index=c(62150,125541,245420,283042,327406),chr=c(1,2,4,5,5))

for (chr in 1:5)

for (chr in 2:5){
png(paste("summary_chr",chr,".png",sep=""))
x<-which(data[,1]==chr)
y.loess<-loess(y~x,span=0.75,data.frame(x=data[x,2],y=lp[x]))
y.predict<-predict(y.loess,data.frame(x=data[x,2]))


#plot(data[x,2],data[x,3],col="grey",type="l",yaxt="n",ylab="")
#axis(4)
#mtext("frequency",side=4,line=3)
#lines(data[x,2],data[x,7],col="black",type="l")
#par(new=T)
plot(data[x,2],lp[x],xlab="pos",ylab="log(p)",main=paste("chr",chr))
lines(data[x,2],y.predict,col="red",type="l")
for(index in qtls[qtls$chr==chr,1]){
 abline(v=data[index,2],col="green")
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


mult_one_fun<-function(locus,r,limit=0.9,cols=c(4,5),memmory=1){
 if((locus%%1000)==0){
  print(locus)
 }
 dist<-abs(shoremap_qtlData[,2]-shoremap_qtlData[locus,2])
 assign("shoremap_q",sapply(dist,q_fun,r),".GlobalEnv")
# assign("shoremap_mod",ppois(dist,r,lower.tail=F),".GlobalEnv")
 assign("shoremap_toUse",shoremap_q>limit,".GlobalEnv")
# assign("shoremap_toUse",ppois(dist,r,lower.tail=F)>0,".GlobalEnv")
 if(sum(shoremap_toUse)>1){
  mle2(ll_one_fun,optimizer="optimize",start=list(f=0.5),fixed=list(locus=locus,r=r),lower=0,upper=1)@coef
 }else{
  shoremap_qtlData[locus,cols[1]]/sum(shoremap_qtlData[locus,cols])
 }
}

ll_one_fun<-function(f,locus,r=50000,cols=c(4,5),memmory=1){
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

