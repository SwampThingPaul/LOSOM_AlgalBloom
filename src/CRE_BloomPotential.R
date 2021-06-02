## 
## Caloosahatchee Algal Bloom risk
##
##
##
## Code was compiled by Paul Julian
## contact info: pjulian@sccf.org

## BAD 
## https://www.tidyverse.org/articles/2017/12/workflow-vs-script/
## Clears Everything...start fresh.
rm(list=ls(all=T));cat("\014");dev.off()

## Libraries
#devtools::install_github("SwampThingPaul/AnalystHelper")
library(AnalystHelper);
library(plyr)
library(reshape2)
library(EnvStats)
library(zoo)

# GIS libraries 
library(rgdal)
library(rgeos)
library(raster)
library(ceramic)

library(flextable)
library(magrittr)

## Paths
wd="C:/Julian_LaCie/_Github/LOSOM_AlgalBloom"

paths=paste0(wd,c("/Plots/","/Export/","/Data/","/src"))
# Folder.Maker(paths);#One and done. Creates folders in working directory.
plot.path=paths[1]
export.path=paths[2]
data.path=paths[3]

GIS.path="C:/Julian_LaCie/_GISData"

# Helper variables
nad83.pro=CRS(SRS_string ="EPSG:4269")
utm17=CRS(SRS_string ="EPSG:26917")
wgs84=CRS(SRS_string = "EPSG:4326")

# GIS ---------------------------------------------------------------------
source("./src/archive/cermanic_key.R")

wmd.mon=spTransform(readOGR(paste0(GIS.path,"/SFWMD_Monitoring_20210119"),"DBHYDRO_SITE_STATION"),utm17)
canals=spTransform(readOGR(paste0(GIS.path,"/SFER_GIS_Geodatabase.gdb"),"SFWMD_Canals"),utm17)
sites=subset(wmd.mon,
       ACTIVITY_S=="Surface Water Grab"&
       STATION%in%c("S79","CES01","S77","PALMOUT","PLN2OUT","TREEOUT"))
plot(sites)

roi=extent(spTransform(gBuffer(sites,width=10000),wgs84))
im <- cc_location(roi,zoom=11)
im2=projectRaster(im,crs=wkt(utm17))
im2=setValues(im2,scales::rescale(values(im2), c(0,255)))

# png(filename=paste0(plot.path,"MonitoringMap.png"),width=6.5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(0.5,0.5,0.5,0.5),oma=c(0.1,0.1,0.1,0.1));
layout(matrix(1:2,2,1,byrow=F),heights=c(1,0.25))
bbox.poly=as(raster::extent(im2),"SpatialPolygons")#makes the polygon
plotRGB(im2)
plot(crop(canals,gBuffer(bbox.poly,width=-750)),add=T,col="lightblue",lwd=1.5,xpd=F)
plot(subset(sites,SITE=="S79"),pch=21,bg="grey",cex=1.25,add=T)
plot(subset(sites,SITE!="S79"),pch=21,bg="indianred1",cex=1.25,add=T)
text(subset(sites,STATION=="CES01"),"STATION",halo=T,hc="black",col="white",pos=4,font=2,cex=0.8)
text(subset(sites,STATION!="CES01"),"STATION",halo=T,hc="black",col="white",pos=2,font=2,cex=0.8)
mapmisc::scaleBar(utm17,"bottomright",outer=T,bty="n",cex=1,seg.len=4,col="white",inset=c(0.01,0.06));

plot(0:1,0:1,ann=F,axes=F,type="n")
legend(0.5,0.5,legend=c("Caloosahatchee Estuary/S79","Lake Okeechobee - Littoral West","Canal"),
       pch=c(21,21,NA),
       lty=c(NA,NA,1),
       lwd=c(0.5,0.5,2),
       pt.bg=c("grey","indianred1",NA),
       col=c("black","black","lightblue"),
       pt.cex=2,ncol=1,cex=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()
# -------------------------------------------------------------------------
dates=date.fun(c("1999-05-01","2021-04-30"))

# Discharge ---------------------------------------------------------------
q.dbkeys=data.frame(SITE=c("S79","S78","S77"),DBKEY=c("00865","DJ236","15635"))
q.cre.dat=data.frame()
for(i in 1:nrow(q.dbkeys)){
  tmp=DBHYDRO_daily(dates[1],dates[2],q.dbkeys$DBKEY[i])
  tmp$DBKEY=as.character(q.dbkeys$DBKEY[i])
  q.cre.dat=rbind(q.cre.dat,tmp)
  print(i)
}
q.cre.dat=merge(q.cre.dat,q.dbkeys,"DBKEY")
q.cre.dat$Date.EST=date.fun(q.cre.dat$Date)

q.cre.dat.xtab=dcast(q.cre.dat,Date.EST~SITE,value.var="Data.Value",mean)
q.cre.dat.xtab$month=as.numeric(format(q.cre.dat.xtab$Date,"%m"))
q.cre.dat.xtab$CY=as.numeric(format(q.cre.dat.xtab$Date,"%Y"))
q.cre.dat.xtab$monCY=with(q.cre.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
q.cre.dat.xtab$WY=WY(q.cre.dat.xtab$Date.EST)


#Positive discharge only
q.cre.dat.xtab$S77=with(q.cre.dat.xtab,ifelse(S77<0,0,S77))
q.cre.dat.xtab$S79=with(q.cre.dat.xtab,ifelse(S79<0,0,S79))
q.cre.dat.xtab$C43=with(q.cre.dat.xtab,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab$Lake=with(q.cre.dat.xtab,ifelse(S77<S79,0,S77-S79))

subset(q.cre.dat.xtab,C43==0&S79>0)

# q.cre.dat.xtab.mon=ddply(q.cre.dat.xtab,c("monCY"),summarise,
#                          S77=sum((S77),na.rm=T),
#                          S79=sum((S79),na.rm=T))
# q.cre.dat.xtab.mon$C43=with(q.cre.dat.xtab.mon,ifelse(S79<S77,0,S79-S77))
# q.cre.dat.xtab.mon$basin.ratio=with(q.cre.dat.xtab.mon,C43/S79)
# # q.cre.dat.xtab.mon$lake.ratio=with(q.cre.dat.xtab.mon,S77/S79)
# q.cre.dat.xtab.mon$Lake=with(q.cre.dat.xtab.mon,ifelse(S77>S79,S79,S77))
# head(q.cre.dat.xtab.mon)


q.cre.dat.xtab.mon=ddply(q.cre.dat.xtab,c("monCY","month",'CY'),summarise,
                         S77=sum(cfs.to.km3d(S77),na.rm=T),
                         S79=sum(cfs.to.km3d(S79),na.rm=T))
q.cre.dat.xtab.mon$C43=with(q.cre.dat.xtab.mon,ifelse(S79<S77,0,S79-S77))
q.cre.dat.xtab.mon$Lake=with(q.cre.dat.xtab.mon,ifelse(S77>S79,S79,S77))
q.cre.dat.xtab.mon$Period=with(q.cre.dat.xtab.mon,ifelse(month%in%c(6:8),"Bloom","Non-Bloom"))

kruskal.test(S77~Period,q.cre.dat.xtab.mon)
kruskal.test(S79~Period,q.cre.dat.xtab.mon)

S79sea.kendall=kendallSeasonalTrendTest(S79~month+CY,data=q.cre.dat.xtab.mon)
print(S79sea.kendall)
S77sea.kendall=kendallSeasonalTrendTest(S77~month+CY,data=q.cre.dat.xtab.mon)
print(S77sea.kendall)

boxplot(S77~month,q.cre.dat.xtab.mon)
boxplot(S79~month,q.cre.dat.xtab.mon)

# png(filename=paste0(plot.path,"BloomRiskPeriod_Qmonth.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,2,0.5,0.5),oma=c(1,2,0.75,0.5));
layout(matrix(1:2,2,1),heights=c(1,0.2))

ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(1,12);xmaj=1:12;xmin=1:12

plot(S77~month,q.cre.dat.xtab.mon,type="n",ylim=ylim.val,xlim=xlim.val,ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
xx=c(5.5,8.5,8.5,5.5)
yy=c(-2,-2,2,2)
polygon(x=xx,y=yy,col=adjustcolor("forestgreen",0.25),border=NA)
with(q.cre.dat.xtab.mon,points(jitter(month-0.2,0.5),S77,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=0.8))
mod=loess(S77~month,q.cre.dat.xtab.mon)
x.val=seq(1,12,length.out=50)
mod.pred=predict(mod,data.frame(month=x.val),se=T)
shaded.range(x.val,mod.pred$fit-qt(0.975,mod.pred$df)*mod.pred$se.fit,mod.pred$fit+qt(0.975,mod.pred$df)*mod.pred$se.fit,bg="grey",lty=0)
lines(x.val,mod.pred$fit)
# lines(x.val,mod.pred$fit-qt(0.975,mod.pred$df)*mod.pred$se.fit,lty=2)
# lines(x.val,mod.pred$fit+qt(0.975,mod.pred$df)*mod.pred$se.fit,lty=2)

with(q.cre.dat.xtab.mon,points(jitter(month+0.2,0.5),S79,pch=21,bg=adjustcolor("indianred1",0.5),col=adjustcolor("black",0.5),lwd=0.01,cex=0.8))
mod=loess(S79~month,q.cre.dat.xtab.mon)
mod.pred=predict(mod,data.frame(month=x.val),se=T)
shaded.range(x.val,mod.pred$fit-qt(0.975,mod.pred$df)*mod.pred$se.fit,mod.pred$fit+qt(0.975,mod.pred$df)*mod.pred$se.fit,bg="indianred1",lty=0)
lines(x.val,mod.pred$fit,col="indianred1")
# lines(x.val,mod.pred$fit-qt(0.975,mod.pred$df)*mod.pred$se.fit,lty=2,col="indianred1")
# lines(x.val,mod.pred$fit+qt(0.975,mod.pred$df)*mod.pred$se.fit,lty=2,col="indianred1")
axis_fun(1,xmaj,xmin,month.abb[1:12],line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=1,line=1.5,"Month")
mtext(side=2,line=2.5,"Discharge (km\u00B3 month\u207B\u00B9)")

plot(0:1,0:1,axes=F,type="n",ylab=NA,xlab=NA);
legend(0.5,-0.5,legend=c("S-77","S79", "LOESS \u00B1 95% CI","Bloom Period"),
       pch=c(21,21,NA,22),lty=c(NA,NA,1,NA),lwd=c(0.1,0.1,1,0),
       col=c(adjustcolor(c("black","black"),0.5),"black",NA),
       pt.bg=c(adjustcolor(c("grey","indianred1"),0.5),NA,adjustcolor("forestgreen",0.25)),
       pt.cex=1.5,ncol=4,cex=0.8,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()


# png(filename=paste0(plot.path,"BloomRiskPeriod.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,2,0.5,0.5),oma=c(2,2,0.75,0.5));
layout(matrix(1:2,1,2))
ylim.val=c(0,0.2);by.y=0.05;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

boxplot(S77~Period,q.cre.dat.xtab.mon,outline=F,ylim=ylim.val,ann=F,axes=F,col=adjustcolor(c("forestgreen","dodgerblue1"),0.5))
axis_fun(1,1:2,1:2,c("High Bloom Risk\n(June - Aug)","Low Bloom Risk\n(Sept - May)"),line=0.25,cex=0.8)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"S-77")
mtext(side=2,line=2.5,"Discharge (km\u00B3 month\u207B\u00B9)")

ylim.val=c(0,1);by.y=0.25;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
boxplot(S79~Period,q.cre.dat.xtab.mon,outline=F,ylim=ylim.val,ann=F,axes=F,col=adjustcolor(c("forestgreen","dodgerblue1"),0.5))
axis_fun(1,1:2,1:2,c("High Bloom Risk\n(June - Aug)","Low Bloom Risk\n(Sept - May)"),line=0.25,cex=0.8)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"S-79")
mtext(side=1,line=1,outer=T,"Bloom Risk Period")
dev.off()

# png(filename=paste0(plot.path,"MonthlyDischarge.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,1.75,0.5,0.5),oma=c(2,2,0.75,0.5));
layout(matrix(1:2,2,1))
ylim.val=c(0,0.6);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=dates;xmaj=seq(xlim.val[1],xlim.val[2],"4 years");xmin=seq(xlim.val[1],xlim.val[2],"1 years")

plot(Lake~monCY,q.cre.dat.xtab.mon,type="n",ylim=ylim.val,xlim=xlim.val,xaxs="i",yaxs="i",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
yr.val=1999:2020
for(i in 1:length(yr.val)){
xx=date.fun(paste0(yr.val[i],c("-06-01","-09-01","-09-01","-06-01")))
yy=c(-10000,-10000,50000,50000)
polygon(x=xx,y=yy,col=adjustcolor("forestgreen",0.25),border=NA)
}
with(q.cre.dat.xtab.mon,shaded.range(monCY,rep(0,length(Lake)),Lake,bg="black",lty=1,lwd=0.25))
# with(q.cre.dat.xtab.mon,lines(monCY,Lake,lwd=1.5))
# with(q.cre.dat.xtab.mon,lines(monCY,C43,lwd=1.5,col="red"))
# axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(1,xmaj,xmin,NA,line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Lake Discharges to CRE")

plot(Lake~monCY,q.cre.dat.xtab.mon,type="n",ylim=ylim.val,xlim=xlim.val,xaxs="i",yaxs="i",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col=adjustcolor("grey",0.5))
for(i in 1:length(yr.val)){
  xx=date.fun(paste0(yr.val[i],c("-06-01","-09-01","-09-01","-06-01")))
  polygon(x=xx,y=yy,col=adjustcolor("forestgreen",0.25),border=NA)
}
with(q.cre.dat.xtab.mon,shaded.range(monCY,rep(0,length(C43)),C43,bg="black",lty=1,lwd=0.25))
axis_fun(1,xmaj,xmin,format(xmaj,"%m-%Y"),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Basin Discharges to CRE")
mtext(side=1,line=2,"Date (Month-Year)")
mtext(side=2,line=0.75,"Discharge (km\u00B3 month\u207B\u00B9)",outer=T)
dev.off()

# WQ Data -----------------------------------------------------------------

params=data.frame(Test.Number=c(18,21,80,20,25,23,61,179,7,16),param=c("NOx","TKN","TN","NH4","TP","SRP","Chla","Chla","Temp","TSS"))
params=subset(params,param%in%c("Temp","Chla"))
wq.sites=data.frame(Station.ID=c("S79","CES01","S77","PALMOUT","PLN2OUT","TREEOUT"),Region=c(rep("CRE",2),rep("LitWest",4)))

wq.dat=DBHYDRO_WQ(dates[1],dates[2],wq.sites$Station.ID,params$Test.Number)
wq.dat=merge(wq.dat,params,"Test.Number")
unique(wq.dat$Collection.Method)
wq.dat=subset(wq.dat,Collection.Method=="G")

wq.dat.xtab=dcast(wq.dat,Station.ID+Date.EST~param,value.var="HalfMDL",mean)
wq.dat.xtab$WY=WY(wq.dat.xtab$Date.EST)
wq.dat.xtab$month=as.numeric(format(wq.dat.xtab$Date.EST,"%m"))
wq.dat.xtab$CY=as.numeric(format(wq.dat.xtab$Date.EST,"%Y"))
wq.dat.xtab$monCY=with(wq.dat.xtab,date.fun(paste(CY,month,"01",sep="-")))
wq.dat.xtab=merge(wq.dat.xtab,wq.sites,"Station.ID")

boxplot(Temp~month,wq.dat.xtab)
abline(h=26)

plot(Chla~Temp,subset(wq.dat.xtab,Region=="CRE"),log='y')
plot(Chla~Temp,subset(wq.dat.xtab,Region=="LitWest"),log='y')

## Regional Comparison
Chla.region=ddply(wq.dat.xtab,c("monCY","month","Region"),summarise,mean.val=mean(Chla,na.rm=T))
Chla.region$bloom=ifelse(Chla.region$mean.val>20,1,0)

# png(filename=paste0(plot.path,"Chl_region.png"),width=6.5,height=4,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1.5,2,0.5,0.5),oma=c(2,2,0.75,0.5));
layout(matrix(1:2,1,2))
ylim.val=c(0,70);by.y=20;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)

boxplot(mean.val~Region,Chla.region,outline=F,ylim=ylim.val,ann=F,axes=F,col=adjustcolor("grey",0.5))
axis_fun(1,1:2,1:2,c("CRE","Littoral West"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,"Chlorophyll-a Concentration (\u03BCg L\u207B\u00B9)")
mtext(side=3,adj=0,"Period of Record - All Months")

boxplot(mean.val~Region,subset(Chla.region,month%in%c(6:8)),outline=F,ylim=ylim.val,ann=F,axes=F,col=adjustcolor("grey",0.5))
axis_fun(1,1:2,1:2,c("CRE","Littoral West"))
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=3,adj=0,"Period of Record - June to August")
mtext(side=1,line=0.5,"Region",outer=T)
dev.off()

kruskal.test(mean.val~Region,Chla.region)
kruskal.test(mean.val~Region,subset(Chla.region,month%in%c(6:8)))
plot(mean.val~month,subset(Chla.region,Region=="CRE"))
plot(mean.val~month,subset(Chla.region,Region=="LitWest"))

## Cross Correlation
Chla.region.xtab=dcast(wq.dat.xtab,monCY~Region,value.var = "Chla",fun.aggregate=function(x) mean(x,na.rm=T))
plot(CRE~LitWest,Chla.region.xtab)

Chla.region.xtab=merge(Chla.region.xtab,q.cre.dat.xtab.mon,"monCY")
# Chla.region.xtab=subset(Chla.region.xtab,as.numeric(format(monCY,"%m"))%in%c(6:8))


layout(matrix(1:4,1,4,byrow=T))
for(h in 0:3){
  tmp.dat=data.frame(CRE.chla=lag(as.zoo(Chla.region.xtab$CRE),0,na.pad=T),
                     Lake.chla=lag(as.zoo(Chla.region.xtab$LitWest),-h,na.pad=T),
                     Lake=Chla.region.xtab$Lake)
  plot(CRE.chla~Lake.chla,tmp.dat,pch=21,bg=ifelse(tmp.dat$Lake>0,"red","grey"))
}

# interesting resources
# https://nwfsc-timeseries.github.io/atsa-labs/sec-tslab-correlation-within-and-among-time-series.html
# https://www.stat.berkeley.edu/~bartlett/courses/153-fall2010/lectures/3.pdf
#
with(subset(Chla.region.xtab,is.na(CRE)==F&is.na(LitWest)==F),ccf(CRE,LitWest,type="correlation"))
# test=with(subset(Chla.region.xtab,is.na(CRE)==F&is.na(LitWest)==F),ccf(CRE,LitWest,type="covariance"))
# abline(h=-1/test$n.used+c(-2,2)/sqrt(test$n.used),col="red")
# abline(h=qnorm((1+0.95)/2)/sqrt(test$n.used))
# x=with(subset(Chla.region.xtab,is.na(CRE)==F&is.na(LitWest)==F),ts.intersect(ts(CRE),ts(LitWest)))
# pacf(x)

ccf.chla.rslt=data.frame()
for(h in 0:15){
  lagged=lag(as.zoo(Chla.region.xtab$LitWest),-h,na.pad=T)
  tmp.dat=as.zoo(Chla.region.xtab$CRE)
  stat=with(data.frame(lag=lagged,dat=tmp.dat),cor.test(lag,dat,method="pearson"))
  ccf.chla.rslt=rbind(ccf.chla.rslt,data.frame(lag=h,estimate=as.numeric(stat$estimate),pval=stat$p.value))
}
ccf.chla.rslt

ylim.val=c(-1,1.1);by.y=0.5;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,15);by.x=5;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/by.x)
# png(filename=paste0(plot.path,"Chla_CCF.png"),width=4,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));
plot(estimate~lag,ccf.chla.rslt,ylim=ylim.val,xlim=xlim.val,type="n",yaxs="i",axes=F,ylab=NA,xlab=NA)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
abline(h=0)
ci.val=qnorm((1+0.95)/2)/sqrt(N.obs(Chla.region.xtab$CRE))
#abline(h=c(ci.val,-ci.val),lty=2,lwd=1.5,col="blue")
polygon(c(-25,25,25,-25),c(ci.val,ci.val,-ci.val,-ci.val),col=adjustcolor("grey",0.5),border=0)
with(ccf.chla.rslt,segments(lag,0,lag,estimate,lwd=1.5,lty=2))
with(ccf.chla.rslt,points(lag,estimate,pch=21,bg=ifelse(pval<0.05,"indianred1","dodgerblue1"),lwd=0.01))
axis_fun(1,line=-0.5,xmaj,xmin,xmaj)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2.5,expression(paste("CCF ",italic("r")["Pearson"])))
mtext(side=1,line=1.5,"Lag (Months)")
legend("topright",legend=c("\u03C1 < 0.05","\u03C1 > 0.05","95% CI"),
       pch=c(21,21,22),pt.bg=c("indianred1","dodgerblue1",adjustcolor("grey",0.5)),col=c("black","black",NA),
       lty=NA,lwd=c(0.1,0.1,0),pt.cex=1.5,cex=0.7,ncol=1,bty="n",y.intersp=1,x.intersp=0.75,xpd=NA,xjust=0.5,yjust=0.5)
dev.off()




bloom.region.xtab=dcast(Chla.region,monCY~Region,value.var = "bloom",fun.aggregate=function(x) mean(x,na.rm=T))
bloom.region.xtab$bothbloom=with(bloom.region.xtab,ifelse(CRE==1&LitWest==1,1,0))
head(bloom.region.xtab)

## 
bloom.region.xtab2=merge(bloom.region.xtab,q.cre.dat.xtab.mon,"monCY")
head(bloom.region.xtab2)

plot(CRE~Lake,bloom.region.xtab2)
plot(CRE~S77,bloom.region.xtab2)
plot(bothbloom~Lake,bloom.region.xtab2)
plot(bothbloom~S77,bloom.region.xtab2)

# library(ggplot2)
# ggplot(bloom.region.xtab2,aes(Lake, CRE))+
#   geom_point(alpha = .15) +
#   geom_smooth(method = "glm", method.args = list(family = "binomial")) +
#   ggtitle("Logistic regression model fit")

mod=glm(CRE~Lake,bloom.region.xtab2,family="binomial")
summary(mod)
confint(mod)
confint.default(mod);# using Standard Errors
# gtsummary::as_flex_table(gtsummary::tbl_regression(mod))

ylim.val=c(0,1);by.y=0.2;ymaj=seq(ylim.val[1],ylim.val[2],by.y);ymin=seq(ylim.val[1],ylim.val[2],by.y/2)
xlim.val=c(0,0.6);by.x=0.2;xmaj=seq(xlim.val[1],xlim.val[2],by.x);xmin=seq(xlim.val[1],xlim.val[2],by.x/2)
# png(filename=paste0(plot.path,"logisticmod_CREChla.png"),width=5,height=3.5,units="in",res=200,type="windows",bg="white")
par(family="serif",mar=c(1,3,0.75,0.75),oma=c(2,1,0.5,0.5));

plot(CRE~Lake,bloom.region.xtab2,type="n",ann=F,axes=F)
abline(h=ymaj,v=xmaj,lty=3,col="grey")
with(bloom.region.xtab2,points(Lake,CRE,pch=21,bg=adjustcolor("grey",0.5),col=adjustcolor("black",0.5),lwd=0.01))

x.val=seq(min(bloom.region.xtab2$Lake,na.rm=T),max(bloom.region.xtab2$Lake,na.rm=T),length.out=200)
mod.pred=predict(mod,data.frame(Lake=x.val), type = "link", se.fit = TRUE)
upr=mod$family$linkinv(mod.pred$fit+(qnorm(0.975)*mod.pred$se.fit))
lwr=mod$family$linkinv(mod.pred$fit-(qnorm(0.975)*mod.pred$se.fit))
fit=mod$family$linkinv(mod.pred$fit)
shaded.range(x.val,lwr,upr,bg="indianred1",lty=0)
lines(x.val,fit,col=adjustcolor("indianred1",0.5),lwd=2)
# lines(x.val,exp(mod.pred[,2]),lty=2)
# lines(x.val,exp(mod.pred[,3]),lty=2)
axis_fun(1,xmaj,xmin,format(xmaj),line=-0.5)
axis_fun(2,ymaj,ymin,format(ymaj));box(lwd=1)
mtext(side=2,line=2,"Algal Bloom (1 = Yes; 0 = No)")
mtext(side=1,line=1.5,"Discharge (km\u00B3 month\u207B\u00B9)")
dev.off()

# diagonistics
plot(1:length(rstandard(mod)),rstandard(mod))
plot(mod,which=4,id.n=5);#identify the top five largest values

# McFadden's r2
mcR2=1-logLik(mod)[1]/logLik(glm(CRE~1,bloom.region.xtab2,family="binomial"))[1]
mcR2
# pscl::pR2(mod)

#classification rate
prop.table(table(bloom.region.xtab2$CRE,predict(mod,data.frame(Lake=bloom.region.xtab2$Lake),type="response")>0.5))
# 83% true positives (based on all data)

data_t <- broom::tidy(mod)
data_g <- broom::glance(mod)
sum_obj <- summary(mod)

data_t$p.value=with(data_t,ifelse(p.value<0.01,"<0.01",ifelse(p.value<0.05,"<0.05",format(round(p.value,2),nsmall=2))))

flextable(data_t, col_keys = c("term", "estimate", "std.error", "statistic", "p.value"))%>%
  colformat_double(j = c("estimate", "std.error", "statistic"), digits = 3)%>%
  # colformat_double(j = c("p.value"), digits = 3)%>%
  set_header_labels(term = "", estimate = "Estimate",
                    std.error = "Standard Error", statistic = "z-value",
                    p.value = "p-value" )%>%
  add_footer_lines(values = c(
    paste("(Dispersion parameter for ", mod$family$family, " family taken to be ", format(sum_obj$dispersion), ")", sep = ""),
    sprintf("Null deviance: %s on %s degrees of freedom", formatC(sum_obj$null.deviance), formatC(sum_obj$df.null)),
    sprintf("Residual deviance: %s on %s degrees of freedom", formatC(sum_obj$deviance), formatC(sum_obj$df.residual)),
    sprintf("McFadden's Pseudo R\u00B2: %s",format(round(mcR2,2),nsmall=2))
    # ,{
    #   if (nzchar(mess <- naprint(sum_obj$na.action)))
    #     paste("  (", mess, ")\n", sep = "")
    #   else character(0)
    # }
  ))%>%
  align(i = 1, align = "right", part = "footer")%>%
  italic(i = 1, italic = TRUE, part = "footer")%>%
  hrule(rule = "auto")%>%
  autofit(part = c("header", "body"))



library(gtsummary)
library(flextable)
library(magrittr)

summary(mod)
as_flex_table(tbl_regression(mod))

#flextable_to_rmd(
flextable::as_flextable(mod)%>%
  compose(i=2,j=1,value=as_paragraph("Lake Discharge"))%>%
  add_footer_lines(value=c("McFadden's R\u00B2 = 0.07"))
#)
# saveRDS(mod,paste0(export.path,"glm_C43.rds"))
# mod=readRDS(paste0(export.path,"glm_C43.rds"))


data.frame(DataType=c(rep("Water Quality",6),rep("Discharge",2)),
           Site=c("S79","CES01","S77","PALMOUT","PLN2OUT","TREEOUT","S77","S79"),
           Region=c(rep("S-79 (CRE)",2),
                    rep("Littoral West (Lake Okeechobee)",4),
                    'From Lake (DBKEY:15635)',"To Estuary (DBKEY: 00865)"))%>%
  flextable()%>%
  merge_v(1)%>%
  fix_border_issues()%>%
  valign(j=1,valign="top")%>%autofit()
