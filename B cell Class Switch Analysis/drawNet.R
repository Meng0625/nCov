setwd("C:/Users/Meng/Desktop/nCov/bCellfNet/")
data = read.csv("./outs/cs_network.csv")
library(plotrix)
data$weight = log2(data$weight+1)

samples_without_bacterial = c("S4","S7","S14","S16", "S12")

samples_with_bacterial = c("S5",  "S6", "S8", "S15", "S9", "S10", "S11", "S1", "S2", "S13")

samples_orderd = c(samples_with_bacterial,samples_without_bacterial)

samples_name_trans = c("P13", "P14",  "P1", "P2", "P6", "P3", "P7",  "P8", "P9",  "P10",  "P16",  "P11",  "P4",  "P12",  "P5")

names(samples_name_trans) = c("S1", "S2",  "S5", "S6", "S7", "S8", "S9", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S4")

Ig.cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")
#[ "IGHM", "IGHG3", "IGHG1", "IGHG2", "IGHA1", "IGHA2"]
gg.list=c('IGHM', 'IGHG3','IGHG1','IGHG2', 'IGHA1', "IGHA2")

# 六个点的横纵坐标
xx.s = c(-2.5,2.5,5,2.5,-2.5,-5)
yy.s = c(4.33,4.33,0,-4.33,-4.33,0)
cc.list = samples_orderd

par(mar=c(0,0,2,0),mfrow=c(3,5))
for(cc in cc.list){
  # sample number transfer
  
  plot(0,0,xlim=c(-6,6),ylim=c(-6,6),axes=F,xlab='',ylab='',cex=0,asp=1,main=samples_name_trans[cc],cex.main=1.5)
  rr.s=c()
  sub_table = subset(data, sample == cc)
  for(ii in 1:length(gg.list)){
    gg1=gg.list[ii] #获得了这个IGHC
    
    # n.gg1=length(which(tmp.CW.mat[,gg1]==1))/10
    # rr.s=c(rr.s,log10(n.gg1))
    for(jj in ii:length(gg.list)){
      if(jj==ii)next
      gg2=gg.list[jj]
      rec =  sub_table[which(sub_table$from==gg1 & sub_table$to ==gg2),]
      if (nrow(rec) != 1){
        print("ERROR")
        print(gg1)
        print(gg2)
        next
      }
      h = rec$weight
      if (rec$qValue< 0.05) {col=rgb(50,50,50,alpha=230,max=255)} else {col=rgb(50,50,50,alpha=40,max=255)}
      segments(xx.s[ii],yy.s[ii],xx.s[jj],yy.s[jj],lwd=h*2,col=col)
    }
  }
  
  for (count in 1:6){
    draw.circle(xx.s[count],yy.s[count],1.3,col = Ig.cols[count],border=rgb(50,50,50,alpha=50,max=255),lwd=3 )
    text(xx.s[count],yy.s[count],sub('IGH','',gg.list[count]),cex=1)
  }
  
  #text(xx.s,yy.s,gg.list,font=2)
}

# 
# plot(0,0,xlim=c(-6,6),ylim=c(-6,6),axes=F,xlab='',ylab='',cex=0,asp=1,cex.main=1.5)
# for (count in 1:3){ 
# draw.circle(x=c(-4,-4,-4)[count],y=c(1,3,5)[count]-0.5,radius=log10(c(20,50,100)/10)[count],col='gray',border=NA)
# text(x=c(-3,-3,-3)[count],y=c(1,3,5)[count]-0.5,c(20,50,100)[count],pos=4,cex=1.5)
# segments(c(-4.5,-4.5,-4.5)[count],c(-5,-3,-1)[count]-0.5, c(-3.5,-3.5,-3.5)[count],c(-5,-3,-1)[count]-0.5,lwd=c(20,50,100)[count]/10,col= rgb(50,50,50,alpha=100,max=255))
# text(x=c(-3,-3,-3)[count],y=c(-5,-3,-1)[count]-0.5,c(20,50,100)[count],pos=4,cex=1.5)
# }


plot(0,0,xlim=c(-6,6),ylim=c(-6,12),axes=F,xlab='',ylab='',cex=0,asp=1,cex.main=1.5)
for (count in 1:6){
  draw.circle(x=-4,y=c(1,3,5,7,9,11)[count]-0.5,radius=1,col=Ig.cols[count],border=NA)
  text(x=-3,y=c(1,3,5,7,9,11)[count]-0.5,gg.list[count],pos=4,cex=1.5)
}
