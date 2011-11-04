#Metod 2... does not work too good


args<-commandArgs(trailingOnly=TRUE)

high_file<-args[1]
low_file<-args[2]
qtl_file<-args[3]
outprefix<-args[4]
ls.t<-read.table(paste(outprefix,".low.AF.txt",sep=""))
hs.t<-read.table(paste(outprefix,".high.AF.txt",sep=""))
#Parameters
minCoverage <- 4
maxCoverage <- 300
chr<-3

#prep data
hs<-read.table(high_file)
ls<-read.table(low_file)

#extract instersection of markers
chrmod<-10**ceiling(log10(max(hs[,2],ls[,2])))
hs.mod<-hs[,1]*chrmod+hs[,2]
ls.mod<-ls[,1]*chrmod+ls[,2]

intersect.mod<-intersect(hs.mod,ls.mod)
hs<-hs[hs.mod %in% intersect.mod,]
ls<-ls[ls.mod %in% intersect.mod,]


hs.cov<-rowSums(hs[,3:5])
ls.cov<-rowSums(ls[,3:5])

goodCov<-(hs.cov>minCoverage & hs.cov<maxCoverage)&(ls.cov>minCoverage &ls.cov<maxCoverage)

hs<-hs[goodCov,]
ls<-ls[goodCov,]

toUse<-hs[,1]==chr
data<-cbind(ls[toUse,2:4],hs[toUse,3:4])



toUse<-ls.t[,1]==chr

#a.t should contain windowed estimates
a.t<-data.frame(pos=ls.t[toUse,2],ls=ls.t[toUse,3]/100,hs=hs.t[toUse,3]/100,)

#use smoothing with zero memmory to get estimates at marker positions
intersect.pos<-intersect(a.t[,1],data[,1])
a.t<-a.t[a.t[,1] %in% intersect.pos,]
data<-data[data[,1] %in% intersect.pos,]

#read qtl file
if(qtl_file==""){
 qtls<-data.frame(pos=-1,chr=-1)
}else{
 qtls<-read.table(qtl_file)
 rank<-sort(sort(abs(qtls[,5]),index.return=TRUE,decreasing=TRUE)$ix,index.return=TRUE)$ix
 lg<-sapply(qtls[,3],function(x) ifelse(sum(qtls[,3]==x)>1,x,0))
 lgType<-sapply(lg,function(x) ifelse(x==0,0,{
  tt<-table(sign(qtls[lg==x,5]))
  ifelse(length(tt)==1,sign(sum(qtls[lg==x,5])),ifelse(tt[1]<tt[2],-tt[1]/tt[2],ifelse(tt[1]==tt[2],0,tt[2]/tt[1])))
 }))
 qtls<-data.frame(id=qtls[,1],pos=qtls[,4],chr=qtls[,3],effect=qtls[,5],rank=rank,lg=lg,lgType=lgType)
}
qtls<-qtls[qtls[,3]==chr,]


freq_sub<-function(x,theta,lambda,effect){
 effect*exp(-4/3*abs(x-theta)/lambda)
}

F=function(x,e1=0.5,e2=-0.5){
 F=matrix(0,nrow=3)
 F[1]= x[3]*x[1]+x[2]-e2
 F[2]= x[3]*x[2]+x[1]-e1
 F[3]= e2*x[1]-e1*x[2]
 F
}

F_wrap<-function(x,e1=0.5,e2=-0.5){
 if(0>=sign(e1)*x[1] || sign(e1)*x[1]>=abs(e1)){
#  print(1)
  NA
 }else if(0>=sign(e2)*x[2] || sign(e2)*x[2]>=abs(e2)){
#  print(2)
  NA
 }else if(abs(x[3]-0.5)>=0.5){
#  print(3)
  NA
 }else{
  sum(F(x,e1,e2)**2)
 }
}


freq<-function(x){
 if(length(x)==1){
  #no QTL linear function between end points
  rep(x[1],length(a.t[,1]))
 }else{
  #edge parameters
  start<-min(a.t[,1])
  end<-max(a.t[,1])+0.01
  b1<-x[1]
  r1<-x[2]
  r2<-x[3]
  #calculate parameters
  qtlest<-matrix(x[-c(1:3)],byrow=TRUE,ncol=2)

  nrOfQtls<-nrow(qtlest)
  par<-c()
  bp<-qtlest[,1]
  if(nrOfQtls>1){
   #sort by position
   qtlest<-qtlest[sort(qtlest[,1],index.return=TRUE)$ix,]
   #initialize
   par<-matrix(c(qtlest[1,1],r1,qtlest[1,2],-1,-1),ncol=5)
   #calculate inbetween parameters
   for(i in 1:(nrOfQtls-1)){
    if(sign(qtlest[i,2])==sign(qtlest[i+1,2])){
     #same direction
     lambda<-4/3*(qtlest[i+1,1]-qtlest[i,1])
     if(qtlest[i,2]>qtlest[i+1,2]){
      #first larger
      lambda<-lambda/log(abs(qtlest[i,2]/qtlest[i+1,2]))
      par<-rbind(par,c(qtlest[i,1],lambda,qtlest[i,2],-1,-1))
     }else{
      lambda<-lambda/log(abs(qtlest[i+1,2]/qtlest[i,2]))
      par<-rbind(par,c(qtlest[i+1,1],lambda,qtlest[i+1,2],-1,-1))
     }
    }else{
     #different directions
     est<-optim(F_wrap,e1=qtlest[i,2],e2=qtlest[i+1,2],par=c(qtlest[i,2]/2,qtlest[i+1,2]/2,0.1),control=list(maxit=10000))$par
     lambda<-4/3*(qtlest[i,1]-qtlest[i+1,1])/log(est[3])
#     x0<-qtlest[i,1]+abs(qtlest[i,2])*(qtlest[i+1,1]-qtlest[i,1])/(abs(qtlest[i,2])+abs(qtlest[i+1,2]))
#     lambda<-4/3*(2*x0-qtlest[i,1]-qtlest[i+1,1])/log(abs(qtlest[i,2]/qtlest[i+1,2]))
#     bp<-c(bp,x0)
#     par<-rbind(par,c(qtlest[i,1],lambda,qtlest[i,2],nrow(par)))
     par<-rbind(par,c(qtlest[i,1],lambda,est[1],qtlest[i+1,1],est[2]))
    }
   }
  }else{
   #par still needs to be initialized
   par<-matrix(c(qtlest[1,1],r1,qtlest[1,2],-1,-1),ncol=5)
  }
  #add last qtl
  par<-rbind(par,c(qtlest[nrOfQtls,1],r2,qtlest[nrOfQtls,2],-1,-1))
  windows<-rowSums(sapply(bp,function(b) ifelse(a.t[,1]<b,0,1)))
  b1+c(sapply(unique(windows),function(win){
   pos<-a.t[windows==win,1]
   p<-par[win+1,]
   if(p[4]<0){
    freq_sub(pos,p[1],p[2],p[3])
   }else{
    freq_sub(pos,p[1],p[2],p[3])+freq_sub(pos,p[4],p[2],p[5])
   }
  }),recursive=TRUE)
 }
}

parameterConstraints<-function(x){
 if(length(x)==1){
  1
 }else if(any(c(x[seq(7,length(x),3)]<min(a.t[,1]),x[seq(7,length(x),3)]>max(a.t[,1])))){
  #position out of range
#  9e5+sum(-x[seq(7,length(x),3)][x[seq(7,length(x),3)]<min(a.t[,1])]+min(a.t[,1]))+sum(x[seq(7,length(x),3)][x[seq(7,length(x),3)]>max(a.t[,1])]-max(a.t[,1]))
  NA
 }else if(any(x[c(2:3,5:6)]<0)){
  #negative edge rates 
#  8e5-sum(x[c(2:3,5:6)][x[c(2:3,5:6)]<0])
  NA
 }else if(any(abs(x[c(seq(8,length(x),3),seq(9,length(x),3))])>1)){
  #too large effects
#  7e5+sum(abs(x[c(seq(8,length(x),3),seq(9,length(x),3))][abs(x[c(seq(8,length(x),3),seq(9,length(x),3))])>1]))
  NA
 }else if(any(abs(x[c(1,4)]-0.5)>0.5)){
  #too high base line
#  6e5+abs(x[1]-0.5)+abs(x[4]-0.5)
  NA
 }else if( any(abs(x[1]+x[c(seq(8,length(x),3))]-0.5)>0.5) ){
  #bad combination of base line and effect
#  5e5+sum(abs(x[c(seq(8,length(x),3))][abs(x[1]+x[c(seq(8,length(x),3))]-0.5)>0.5]))
  NA
 }else if( any(abs(x[4]+x[c(seq(9,length(x),3))]-0.5)>0.5) ){
  #bad combination of base line and effect
#  5e5+sum(abs(x[c(seq(9,length(x),3))][abs(x[4]+x[c(seq(9,length(x),3))]-0.5)>0.5]))
  NA
 }else{
  1
 }
}

estLambda<-function(x,col=2){
 f.est<-freq(x)
 sum((a.t[,col]-f.est)**2)
}

separateParameters<-function(x){
 m<-matrix(x,nrow=3)
 base<-m[,1:2]
 qtlest<-matrix(m[,-(1:2)],nrow=3)
 x1<-c(base[,1],qtlest[1:2,])
 x2<-c(base[,2],qtlest[c(1,3),])
 data.frame(x1=x1,x2=x2)
}

freqBoth<-function(x){
 pc<-parameterConstraints(x)
 if(is.na(pc)){
  pc
 }else{
  if(length(x)==1){
   estLambda(x,2)+estLambda(x,3)
  }else{
   par<-separateParameters(x)
   estLambda(par$x1,2)+estLambda(par$x2,3)
  }
 }
}

llBoth<-function(x){
 pc<-parameterConstraints(x)
 if(is.na(pc)){
  pc
 }else{
  f1<-c()
  f2<-c()
  if(length(x)==1){
   f1<-freq(x)
   f2<-freq(x)
  }else{
   par<-separateParameters(x)
   f1<-freq(par$x1)
   f2<-freq(par$x2)
  }
  d1<-cbind(data[,2],rowSums(data[,2:3]),f1)
  d2<-cbind(data[,4],rowSums(data[,4:5]),f2)
  -sum(apply(d1,1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))-sum(apply(d2,1,function(r) dbinom(r[1],size=r[2],prob=r[3],log=TRUE)))
 }
}

calcbic<-function(x){
 2*llBoth(x)+length(x)*log(nrow(data))
}

estPlot<-function(x){
 plot(a.t[,c(1,3)],type="l",ylim=c(0,1))
 lines(a.t[,1:2])
 apply(matrix(unique(qtls[,2]),ncol=1),1,function(x) abline(v=x))
 par<-separateParameters(x)
 lines(a.t[,1],freq(par$x1),col="red")
 lines(a.t[,1],freq(par$x2),col="red")
 apply(matrix(unique(x[seq(7,length(x),3)]),ncol=1),1,function(x) abline(v=x,col="red"))
}

ll1w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9)
 l<-llBoth(x)
 print(l)
# l<-round(l*1000)/1000
 l
}

ll2w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12)
 l<-llBoth(x)
# print(l)
# l<-round(l*1000)/1000
 l
}

ll3w<-function(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15){
 x<-c(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15)
 l<-llBoth(x)
# print(l)
# l<-round(l*1000)/1000
 l
}



par.start<-c(0.4,20e6,20e6,0.6,20e6,20e6,3e7,-0.2,0.2)

system.time({
opt<-optim(freqBoth, par=par.start, method="Nelder-Mead", control=list(parscale=par.start, maxit=500))
print(opt$val)
val<-Inf
while (opt$val<val-1e-4 && opt$convergence==1) {
 val<-opt$val
 opt<-optim(freqBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=opt$par, maxit=500))
 print(opt$val)
}
})

par.start<-opt$par

system.time({
opt<-optim(llBoth,par=par.start, method="Nelder-Mead", control=list(parscale=par.start))
print(opt$val)
val<-Inf
while (opt$val<val && opt$convergence==1) {
#while (opt$convergence==1) {
 val<-opt$val
 opt<-optim(llBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=opt$par))
 print(opt$val)
}
})

fixed<-c(list(),ci@coef)

mle.start<-c(list(),m@coef)

mle.start<-c(list(),opt$par)
names(mle.start)<-paste("p",1:length(opt$par),sep="")

system.time(m<-mle2(ll1w, start=mle.start))
system.time(ci<-confint(m,parm=c(7)))

system.time(ci<-confint(m,parm=c(7),control=list(parscale=opt$par)))






















# to reimplement
calculateScale<-function(x){
 scaleBaseline<-1
 scaleQTL<-c(1,1,1,1,1,1,1)
 (x/1)*c(rep(scaleBaseline,2),rep(scaleQTL,(length(x)-2)/length(scaleQTL)))
}

addQTL<-function(x,direction=1){
 direction=sign(direction)
 if(direction==0){
  direction<-1
 }
 par.matrix.both<-parameterToMatrix(x)
 par.matrix.1<-matrix(par.matrix.both[,1:5], ncol=5)
 par.matrix.2<-matrix(par.matrix.both[,6:10], ncol=5)
 baseLine<-x[c(1:2)]
 curEst<-x[-c(1:2)]
 nrOfQtls<-length(curEst)/7
 #adjust the baseline
 baseLine<-baseLine*nrOfQtls/(nrOfQtls+1)
 curMat<-matrix(curEst,byrow=TRUE,ncol=7) 
 #position the new QTL in the middle of the gap with the highest average deviation from the frequency
 mn<-min(a.t[,1])
 mx<-max(a.t[,1])
 d<-diff(c(mn,curMat[,7],mx))
 f1<-freq_tot(par.matrix.1)
 f2<-freq_tot(par.matrix.2)
 section<-rowSums(sapply(curMat[,7],function(x) ifelse(a.t[,1]<x,0,1)))
 deviation<-tapply(1:length(section),section,function(i) mean(abs(a.t[i,2]-f1[i])**2+abs(a.t[i,3]-f2[i])**2))
 pos<-mn+sum(d[1:which.max(deviation)])-d[which.max(deviation)]/2
 #new QTL
 q<-c(apply(matrix(curMat[,1:6],ncol=6),2,mean),pos)
 q[5:6]<-direction*q[5:6]
 #use the mean of the previous qtls for the other parameters
 c(baseLine,curEst,q)
}



method<-"Nelder-Mead"
qtl_number<-1
direction1<--sign(diff(mean(a.t[,2:3])))/2
qtl_start <-c(rep(max(a.t[,1])/4/qtl_number,4),direction1/qtl_number,-direction1/qtl_number,max(a.t[,1])/(qtl_number+1))
par.start<-c(rep(0.5/qtl_number,2),rep(qtl_start,qtl_number)*c(sapply(1:qtl_number, function(x) c(rep(1,6),x))))
#scale<-calculateScale(par.start)
#par.start<-par.start/scale




system.time({
opt<-optim(l1Est, par=par.start, method="Nelder-Mead", control=list(parscale=calculateScale(par.start), maxit=1000))
print(opt$val)
val<-Inf
while (opt$val<val-1e-4 && opt$convergence==1) {
 val<-opt$val
 opt<-optim(l1Est, par=opt$par, method="Nelder-Mead", control=list(parscale=calculateScale(opt$par), maxit=1000))
 print(opt$val)
}
})

par.start<-opt$par

system.time({
opt<-optim(llBoth,par=par.start, method="Nelder-Mead", control=list(parscale=calculateScale(par.start)))
print(opt$val)
val<-Inf
#while (opt$val<val && opt$convergence==1) {
while (opt$convergence==1) {
 val<-opt$val
 opt<-optim(llBoth, par=opt$par, method="Nelder-Mead", control=list(parscale=calculateScale(opt$par)))
 print(opt$val)
}
})

mle.start<-c(list(),ci@coef)

mle.start<-c(list(),opt$par)
names(mle.start)<-paste("p",1:length(opt$par),sep="")

system.time(m<-mle2(ll2w, start=mle.start,control=list(parscale=calculateScale(opt$par),reltol=1e-3)))
system.time(ci<-confint(m,parm=c(9,16),control=list(parscale=calculateScale(opt$par),reltol=1e-3)))









