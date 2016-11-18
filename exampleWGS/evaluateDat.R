#This is to investigate what kind of filtering should be done with swainyPrepperWGS.R
library("ggplot2")
library(dplyr)

setwd("/home/mac/swainysmoother/exampleWGS/test")

#mac@swc-megabox-dt:~/swainysmoother/exampleWGS/test$ grep omy01 M027288.dat > M027288.omy01.dat
#mac@swc-megabox-dt:~/swainysmoother/exampleWGS/test$ grep omy01 M027304.dat > M027304.omy01.dat

data1<-read.table("M027288.omy01.dat")
ranges1<-data1$V3-data1$V2
data1$ranges<-data1$V3-data1$V2

data2<-read.table("M027304.omy01.dat")
ranges2<-data2$V3-data2$V2
data2$ranges<-data2$V3-data2$V2

data3<-read.table("M075289.omy01.dat")
ranges3<-data3$V3-data3$V2
data3$ranges<-data3$V3-data3$V2

data4<-read.table("M075304.omy01.dat")
ranges4<-data4$V3-data4$V2
data4$ranges<-data4$V3-data4$V2


mean(ranges1)
mean(ranges2)
mean(ranges3)
mean(ranges4)
max(ranges1)
max(ranges2)
max(ranges3)
max(ranges4)

#df1<-as.data.frame(data1)
#df2<-as.data.frame(data2)
#df1$join<-'ind1'
#df2$join<-'ind2'
#df12<-rbind(df1,df2)
#ggplot(df12, aes(ranges, fill=join))+geom_density(alpha=0.3)

df1<-as.data.frame(subset(data1, ranges<500))
df2<-as.data.frame(subset(data2, ranges<500))
df3<-as.data.frame(subset(data3, ranges<500))
df4<-as.data.frame(subset(data4, ranges<500))
df1$join<-'ind1'
df2$join<-'ind2'
df3$join<-'ind3'
df4$join<-'ind4'
df1234<-rbind(df1,df2,df3,df4)

pdf("./coveragePlot.pdf")
ggplot(df1234, aes(ranges, fill=join))+geom_density(alpha=0.3)
dev.off()
