# Read in command line options
args <- commandArgs()

pdf(file=args[5], width=30, height=10)
#layoutmat=matrix(data=c(1,2,3,4,5), ncol=1, nrow=5)
#layout(layoutmat)

options(scipen=999999999)

data<-read.table(args[6])


# Plot
plot(1,1, xlim=c(min(data$V2), max(data$V2)), ylim=c((-1)*max(data$V4), max(data$V3)), type="n", axes=F, xlab=data$V1[1], ylab="Read count")
segments(data$V2, 0, data$V2, data$V3, col="red")
segments(data$V2, 0, data$V2, data$V4*(-1), col="blue")
axis(2)
axis(1)



# Thank you for riding:
dev.off()

