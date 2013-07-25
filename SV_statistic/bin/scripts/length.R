pdf("HEG4.SV.length.pdf")

barlegend <- function(x,y,height,length,name,color){
    rect(x,y-0.04*height,x+0.05*length,y,col=color,border=FALSE)
    text(x+0.15*length,y-0.02*height,labels=name)
}


par(mar=c(6,4,4,2))
x <- read.table("HEG4.SV.length.distr",skip=1)
data <- rbind(x[,3]/sum(x[,3]),x[,4]/sum(x[,4]))
xx <- barplot(data,beside=TRUE,ylab="Proportion",border=FALSE,ylim=c(0,0.7),col=c("Orange","blue"))
axis(1,c(0.5,max(xx)+0.5),line=0,labels=c("",""))
text(xx[1,]+0.2,rep(-0.08,6),offset=2,labels=x[,2],srt=55,xpd=TRUE)
text(xx[1,],x[,3]/sum(x[,4])+0.03,offset=2,labels=x[,3],srt=55,xpd=TRUE)
text(xx[1,]+1.2,x[,4]/sum(x[,4])+0.03,offset=2,labels=x[,4],srt=55,xpd=TRUE)
legend("topright",c("Deletion","Insertion"),bty="n",border="NA",lty=c(0,0),cex=1,fill=c("orange","blue"))
#barlegend(11,0.47,0.6,16,"Deletion","Orange")
#barlegend(11,0.44,0.6,16,"Insertion","Blue")
mtext("Indel Length (bp)",side=1, at=13,line=5)

dev.off()

