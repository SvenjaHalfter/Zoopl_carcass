#Analysis for carcass publication

#set working directory
#setwd("C:/Users/shalfter/Documents/Outputs/2 Under_review/Carcasses/R/For_Github")

#open libraries
library(ncdf4) #opening ncdf files
library(tidyverse) #plotting and data wrangling
library(reshape2) #melt function
library(sp)#working with spatial data (Spatialpointsdataframe function)
library(raster)#working with spatial data (raster function)
library(mapdata)#get world data
library(maps)#required for mapdata
library(maptools) #working with spatial data (map2SpatialPolygons function)
library(rasterVis)#working with spatial data (rasterTheme, levelplot function)
library(viridis)#colour scheme for map
library(cowplot)#plot grids
library(lme4) #linear model

####Figure 1 - Map####
#code by Emma Cavan (https://github.com/e-cavan/Southern_Ocean_Microbial_Respiration/blob/master/Cavan_Boyd_2018.R)
#downloaded the satellite chlorophyll file from https://oceancolor.gsfc.nasa.gov/cgi/l3 - September/October 2018
#select Aqua MODIS Chl at 4km resolution

#aim: extract ocean colour from netcdf file, transfer into csv and plot in ggplot2
#study area 
lonmax<-155 
lonmin<-140 
latmax<--40 
latmin<--50 

file.exists("MODISA_chl.csv")#check if a csv file exists in working directory:
# FALSE: continue below
# TRUE: skip the data extraction and jump directly to plotting
# Create list (not necessary here, but useful when you have several files)
f <- list.files(".", pattern="*.L3m_MO_CHL_chlor_a_4km.nc",full.names=F)
var<-"chlor_a"# variable of interest
# open netCDF file              
data<-nc_open(f[1])
# extract data
lon<-ncvar_get(data,"lon")
lat<-ncvar_get(data,"lat")
value<-ncvar_get(data,var)
unit<-ncatt_get(data,var,"units")$value
# matrix to data.frame
dimnames(value)<-list(lon=lon,lat=lat)
dat.var<-melt(value,id="lon")
# select data from the study area taking out missing data
dat.varSAtmp<-subset(dat.var,lon<=lonmax & lon>=lonmin & lat<=latmax &
                       lat>=latmin & value<45)
# extract date information
dateini<-ncatt_get(data,0,"time_coverage_start")$value
dateend<-ncatt_get(data,0,"time_coverage_end")$value
datemean<-mean(c(as.Date(dateend,"%Y-%m-%dT%H:%M:%OSZ"),as.Date(dateini,"%Y-%m-%dT%H:%M:%OSZ")))
year<-substring(datemean,0,4)
month<-substring(datemean,6,7)
# prepare final data set
dat.varSA<-data.frame(rep(as.integer(year,nrow(dat.varSAtmp))),rep(as.integer(month,nrow(dat.varSAtmp))),
                      dat.varSAtmp,rep(unit,nrow(dat.varSAtmp)),rep(var,nrow(dat.varSAtmp)))
names(dat.varSA)<-c("year","month","lon","lat","value","unit","var")
# save csv file
fe<-file.exists("MODISA_chl.csv")
write.table(dat.varSA,"MODISA_chl.csv",row.names=FALSE,col.names=!fe,sep=",",dec=".",append=fe)
# close connection
nc_close(data)

#### Start here if csv already exists
# Convert dataframe to raster for plotting
chl <- read.csv('MODISA_chl.csv', sep=',')
chl_m<-chl[,c(3,4,5)]# subset data - lon, lat, value of 1 plot
chl_m <- subset(chl_m, select=c(3,1:2))# put value (chl) column first
chl_m <- chl_m[chl_m$value<2,]# set max chl level
xy<-chl_m[,c(2,3)]# data frame of lon and lat
# convert to spatial data frame fram using value df and lon/lat df
spdf <- SpatialPointsDataFrame(coords = xy, data = chl_m,
                               proj4string = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
# convert pixels due to irregular point data, tolerance given from gridden(pixels)=TRUE
pixels <- SpatialPixelsDataFrame(spdf, tolerance = 0.00448656, spdf@data)
gridded(pixels) = TRUE# tell r it is gridded
# convert to raster using raster package
r<-raster(pixels)

# levelplot 
boundaries<-map("worldHires", ylim=c(-50,-40), xlim=c(140,155), col='gray90', 
                fill=TRUE, myborder=0.01, boundary=FALSE)
IDs <- sapply(strsplit(boundaries$names, ":"), function(x) x[1])
bPols <- map2SpatialPolygons(boundaries, IDs=IDs,
                             proj4string=CRS(projection(r)))
mapTheme <- rasterTheme(region=viridis(300))
b <- levelplot(r, contour=FALSE, labels=FALSE, margin=FALSE, par.settings=mapTheme,
               ylim=c(-50,-40), xlim=c(140,155),family="A", 
               scales=list(x=list(cex=2),y=list(cex=2)), xlab=list("Longitude", cex=2),
               ylab=list("Latitude", cex=2),
               main=list(expression(paste('Chl a (',mu,'g L'^'-1'*')')),cex=2),
                         colorkey=list(at=seq(0, 2, 0.02), labels=list(at=c(0,1,2), 
                                                             labels=c('0','1','2'), 
                                                             cex.lab=2)))
par(mar=c(1,1,1,1))
# coordinates of Process station 1 and 2
t <- c(141.5,153) #lon
y<- c(-46,-45.5) #lat
xy <- data.frame(t,y)
yx <- SpatialPoints(xy[,1:2])
yx$ID<-c("PS1", "PS2")
# plot map with Tasmania and sampling sites (orange points)
mapSOTS <-b + layer(sp.polygons(bPols, fill='darkslategray')) + 
  layer(sp.points(yx, pch=16, col='darkorange2', cex=2), columns = 1)+
  layer(sp.text(coordinates(yx),txt=yx$ID, col="white",pos=1, cex=2))
mapSOTS

####Figure 1 - CTD####

poc<- read.csv('POM_measured.csv', header=T) #POC/PON data from Niskin 
ctd1 <- read.csv('CTD_stationPS1.csv', header=T) #CTD data station PS1
ctd2 <- read.csv('CTD_stationPS2.csv', header=T) #CTD data station PS2
chl <- read.csv('Chlorophyll_measured.csv', header=T) #Chlorophyll data from Niskin 


#Plotting PS1
#Using base R instead of ggplot2 because it deals better with several axes
par(mfrow=c(1,3), mar=c(12,11,7,1)) #extra margin room mar(bottom, left, top, right)

plot(poc$CperL[poc$Station=="PI"], poc$Depth[poc$Station=="PI"],type= 'b', pch=16,
     col="darkslategray", lwd=2, cex=2,ylim=c(200,0), ylab="", xlab='', 
     xaxt='n', yaxt='n', cex.axis=2, cex.lab=2) 
title(ylab='Depth (m)', mgp=c(4,1,0),cex.lab=2)
axis(3, at=seq(0,30,5), cex.axis=2)
axis(2, labels=FALSE)
mtext(expression(paste('POC (',mu,'g L'^'-1'*')')), 3, line = 3.5, cex=1.5,
      col="darkslategray")
text(y=seq(200, 0, -50),cex=2,
     par('usr')[1], pos=2,labels=as.vector(c(200 , 150 , 100 , 50 ,0 )),srt=0, 
     xpd=TRUE)
# add PON
par(new=T)
plot(poc$NperL[poc$Station=="PI"], poc$Depth[poc$Station=="PI"], type='b', pch=17,
     col="darkslategray4",lwd=2, cex=2,axes=F, xlab='', ylab='', ylim=c(200,0), 
     xlim=c(2.5,5))
axis(1, at=seq(2.5,5,0.5), cex.axis=2, col="darkslategray4", col.axis="darkslategray4")
mtext(expression(paste('PON (',mu,'g L'^'-1'*')')), 1, line=3.5, cex=1.5, 
      col="darkslategray4")
# add chlorophyll
par(new=T)
plot(chl$Chl..ug.l.[chl$CTD.ID=='CTD 23'], chl$Depth[chl$CTD.ID=='CTD 23'], 
     type='b', col="coral", pch=18, lwd=2,cex=2, axes=F, xlab='', ylab='', ylim=c(200,0),
     xlim=c(0, 0.5))
axis(1, at=seq(0,0.5, 0.05), cex.axis=2, line=5, col="darkorange2", 
     col.axis="darkorange2")
mtext(expression(paste('Chlorophyll (',mu,'g L'^'-1'*')')), 1, line=9, cex=1.5, 
      col="darkorange2")
#add euphotic zone
abline(h=126, col="black", lty=2)

#Plotting PS2
plot(poc$CperL[poc$Station=="PII"], poc$Depth[poc$Station=="PII"],type= 'b', pch=16,
     col="darkslategray", lwd=2, cex=2,ylim=c(200,0), ylab="", xlab='', 
     xaxt='n', yaxt='n', cex.axis=2, cex.lab=2) 
title(ylab='Depth (m)', mgp=c(4,1,0),cex.lab=2)
axis(3, at=seq(0,30,5), cex.axis=2)
axis(2, labels=FALSE)
mtext(expression(paste('POC (',mu,'g L'^'-1'*')')), 3, line = 3.5, cex=1.5,
      col="darkslategray")
text(y=seq(200, 0, -50),cex=2,
     par('usr')[1], pos=2,labels=as.vector(c(200 , 150 , 100 , 50 ,0 )),srt=0, 
     xpd=TRUE)
# add PON
par(new=T)
plot(poc$NperL[poc$Station=="PII"], poc$Depth[poc$Station=="PII"], type='b', pch=17,
     col="darkslategray4",lwd=2, cex=2,axes=F, xlab='', ylab='', ylim=c(200,0), 
     xlim=c(2.5,5))
axis(1, at=seq(2.5,5,0.5), cex.axis=2, col="darkslategray4", col.axis="darkslategray4")
mtext(expression(paste('PON (',mu,'g L'^'-1'*')')), 1, line=3.5, cex=1.5, 
      col="darkslategray4")
# add chlorophyll
par(new=T)
plot(chl$Chl..ug.l.[chl$CTD.ID=='CTD 19'], chl$Depth[chl$CTD.ID=='CTD 19'], 
     type='b', col="coral", pch=18, lwd=2,cex=2, axes=F, xlab='', ylab='', ylim=c(200,0),
     xlim=c(0, 0.5))
axis(1, at=seq(0,0.5, 0.05), cex.axis=2, line=5, col="darkorange2", 
     col.axis="darkorange2")
mtext(expression(paste('Chlorophyll (',mu,'g L'^'-1'*')')), 1, line=9, cex=1.5, 
      col="darkorange2")
#add euphotic zone
abline(h=116, col="black", lty=2)


####Environmental conditions 0-200 m####
#Summarise temperature and salinity
ctd1%>%
  filter(Pressure<=200)%>%
  summarize(Temp=mean(Temperature, na.rm=T), Temp.SD=sd(Temperature,na.rm=T),
            Sal=mean(Salinity, na.rm=T),Sal.SD=sd(Salinity,na.rm=T))
ctd2%>%
  filter(Pressure<=200)%>%
  summarize(Temp=mean(Temperature, na.rm=T), Temp.SD=sd(Temperature,na.rm=T),
            Sal=mean(Salinity, na.rm=T),Sal.SD=sd(Salinity,na.rm=T))
#calculate depth of EZ (1 % PAR)
ctd1%>%
  filter(Biospherical.PAR<=1.56)
ctd2%>%
  filter(Biospherical.PAR<=0.123)

#calculate MLD (following Hosoda 2010, delta T=0.2 C)
ctd1%>%
  filter(Temperature<=9.227)
ctd2%>%
  filter(Temperature<=10.415)

#calculate mean POC/PON in the upper 200 m
poc%>%
  filter(Station=="PI")%>%
  summarize(POC_mean=mean(CperL, na.rm=T), POC_sd=sd(CperL, na.rm=T),
            PON_mean=mean(NperL, na.rm=T), PON_sd=sd(NperL, na.rm=T))
poc%>%
  filter(Station=="PII")%>%
  summarize(POC_mean=mean(CperL, na.rm=T), POC_sd=sd(CperL, na.rm=T),
            PON_mean=mean(NperL, na.rm=T), PON_sd=sd(NperL, na.rm=T))

####Figure 2 - N.tonsus size distribution####
length <- read.csv('Prosome_length_sample.csv', header=T)

Species1<-ggplot(length)+
  geom_histogram(aes(x=PS1, fill=Daytime), position="dodge", show.legend=FALSE)+
  scale_fill_manual(name= "Sampling", labels= c("Day", "Night"), values=c("coral", "darkslategray"))+
  labs(title="PS1", x='',  y="Frequency")+
  geom_vline(xintercept=1.93, linetype="dashed", size=1)+
  geom_vline(xintercept=2.61, linetype="dashed", size=1)+
  annotate(geom="text", x=1.6, y=60, label="CIV", size=8)+
  annotate(geom="text", x=2.25, y=60, label="CV", size=8)+
  annotate(geom="text", x=3.3, y=60, label="CVI (adult)", size=8)+
  theme_cowplot(12)+
  theme(text=element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size=20))+
  coord_cartesian(xlim=c(1.5,4.0), ylim=c(0,60))

Species2<-ggplot(length)+
  geom_histogram(aes(x=PS2, fill=Daytime), position="dodge", show.legend=FALSE)+
  scale_fill_manual(name= "Sampling", labels= c("Day", "Night"), values=c("coral", "darkslategray"))+
  labs(title="PS2", x='Prosome length (mm)',  y="Frequency")+
  geom_vline(xintercept=1.93, linetype="dashed", size=1)+
  geom_vline(xintercept=2.61, linetype="dashed", size=1)+
  annotate(geom="text", x=1.6, y=60, label="CIV", size=8)+
  annotate(geom="text", x=2.25, y=60, label="CV", size=8)+
  annotate(geom="text", x=3.3, y=60, label="CVI (adult)", size=8)+
  theme_cowplot(12)+
  theme(text=element_text(size = 20), axis.text = element_text(size=20), 
        axis.title = element_text(size=20))+
  coord_cartesian(xlim=c(1.5,4.0), ylim=c(0,60))

legend.plot <- ggplot(length)+
  geom_histogram(aes(x=PS1, fill=Daytime), position="dodge")+
  scale_fill_manual(name= "Sampling", labels= c("Day", "Night"), values=c("coral", "darkslategray"))+
  labs(title="PS1", x='',  y="Frequency")+
  theme_cowplot(12)+
  theme(text=element_text(face="bold", size=22))+
  coord_cartesian(xlim=c(1.5,4.0), ylim=c(0,60))

right_side <- plot_grid(Species1, Species2, labels=c('a', 'b'), label_size=20,ncol=1)

legend <- get_legend(
  legend.plot + theme(legend.box.margin = margin(0, 0, 0, 14)))
plot_grid(right_side, legend, rel_widths = c(3, .4))

####Carbon standing stock in N.tonsus####
#Using the equations calculated in figure 3

length$PS1_carbon<-(41.638*length$PS1)+23.523 #carbon based on prosome length
length$PS2_carbon<-(41.638*length$PS2)+23.523
Cstock<-length%>%
  group_by(Daytime)%>%
  summarize(PS1= sum(PS1_carbon,na.rm=T),PS2=sum(PS2_carbon,na.rm=T)) 
#result: total C stock day/night in ug
Cstock$PS1<-(Cstock$PS1/1000)*2  #ug to mg and accounting for 0.5 m2 diameter of the net
Cstock$PS2<-(Cstock$PS2/1000)*2




####Figure 3 - CHN & sinking speed plots####
carbon <- read.csv("Copepod_carbon.csv", header = T)
speed<- read.csv("Sinking_speed.csv", header = T)

#A: dry weight versus carbon
carbon$dry.weight.ug<-carbon$dry.weight..mg.*1000

c<-lm((carbon$C..ug.)~(carbon$dry.weight.ug)); summary(c) 
carbon_mod<-seq(0,0.6,0.1)
met<-(0.115*carbon_mod)+112.930 #put in calculated m and b

plotA<-ggplot(data=carbon, aes(dry.weight.ug, C..ug.))+
  geom_point(size=2)+
  lims(y=c(90,200))+
  labs(x=expression(paste('Dry weight (',mu,'g)')), y=expression(paste('Carbon (',mu,'g)')))+
  geom_smooth(method="lm", se=F, colour="darkslategray4")+
  annotate(geom="text", label="p < 0.01", x=500,y=90, size=6)+
  annotate(geom="text", label=expression(paste('R'^'2'=='0.69')), x=500, y=100, size=6)+
  theme_cowplot(12)+
  theme(text=element_text(size=16), axis.text = element_text(size=20),
        axis.title = element_text(size=20))

#B: prosome length versus carbon 
c<-lm((carbon$C..ug.)~(carbon$length..mm.)); summary(c) 
carbon_mod<-seq(0,0.6,0.1)
met<-41.638*carbon_mod+23.523 #put in calculated m and b

plotB<-ggplot(data=carbon, aes(length..mm., C..ug.))+
  geom_point(size=2)+
  lims(x=c(2.8,3.6),y=c(90,200))+
  labs(x="Prosome length (mm)", y=expression(paste('Carbon (',mu,'g)')))+
  geom_smooth(method="lm", se=F, colour="darkslategray4")+
  annotate(geom="text", label="p < 0.01", x=3.45,y=90, size=6)+
  annotate(geom="text", label=expression(paste('R'^'2'=='0.24')), x=3.45, y=100, size=6)+
  theme_cowplot(12)+
  theme(text=element_text(size=16), axis.text = element_text(size=20),
        axis.title = element_text(size=20))

#C: sinking speed 
r<-lm((speed$speed..m.min.)~(speed$size..mm.)); summary(r)
size_mod<-seq(1,4,0.5)
met<-(0.21454*size_mod)-0.18047 #put in calculated m and b

plotC<-ggplot(data=speed, aes(size..mm., speed..m.min.))+
  geom_point(size=2)+
  lims(x=c(2.8,3.6),y=c(0,0.8))+
  labs(x="Prosome length (mm)", y=expression(paste('Sinking velocity (m min '^'-1'*')')))+
  geom_smooth(method="lm", se=F, colour="darkslategray4")+
  annotate(geom="text", label="p < 0.01", x=3.45,y=0.0, size=6)+
  annotate(geom="text", label=expression(paste('R'^'2'=='0.07')), x=3.45, y=0.075, size=6)+
  theme_cowplot(12)+
  theme(text=element_text(size=16), axis.text = element_text(size=20),
        axis.title = element_text(size=20))


plot_grid(plotA, plotB, plotC, labels= 'auto', nrow=1, label_size = 20)

####Figure 4 - Decomposition####
# Code by Emma Cavan 
#experiment on the 29092018 - batch 1
setwd("~/Outputs/2 Under_review/Carcasses/R/For_Github/batch_1")
t0 <- read.csv("t0.csv", header=T)
t1 <- read.csv("t1.csv", header=T)
t2 <- read.csv("t2.csv", header=T)
t3 <- read.csv("t3.csv", header=T)
t4 <- read.csv("t4.csv", header=T)
t5 <- read.csv("t5.csv", header=T)
t6 <- read.csv("t6.csv", header=T)
t7 <- read.csv("t7.csv", header=T)
t8 <- read.csv("t8.csv", header=T)

# create one data frame with all rows - all dates and stations
resp1<-rbind(t0, t1, t2, t3, t4, t5, t6, t7, t8)
# Add chamber column 
A01<-rep('A01', 4)
A02<-rep('A02', 4)
A03<-rep('A03', 4)
A04<-rep('A04', 4)
B05<-rep('B05', 4)
B06<-rep('B06', 4)
B07<-rep('B07', 4)
B08<-rep('B08', 4)
C09<-rep('C09', 4)
C10<-rep('C10', 4)
C11<-rep('C11', 4)
C12<-rep('C12', 4)
chamber<-c(A01, A02, A03, A04, B05, B06, B07, B08, C09, C10, C11, C12)
resp1$chamber<-rep(chamber)
# Add t column
T0 <- rep('T0', 48)
T1 <- rep('T1', 48)
T2 <- rep('T2', 48)
T3 <- rep('T3', 48)
T4 <- rep('T4', 48)
T5 <- rep('T5', 48)
T6 <- rep('T6', 48)
T7 <- rep('T7', 48)
T8 <- rep('T8', 48)
resp1$t <- c(T0, T1, T2, T3, T4, T5, T6, T7,T8)

#experiment on the 01/10/2018 - batch 2
setwd("~/Outputs/2 Under_review/Carcasses/R/For_Github/batch_2")
t0 <- read.csv("t0.csv", header=T)
t1 <- read.csv("t1.csv", header=T)
t2 <- read.csv("t2.csv", header=T)
t3 <- read.csv("t3.csv", header=T)
t4 <- read.csv("t4.csv", header=T)
t5 <- read.csv("t5.csv", header=T)
t6 <- read.csv("t6.csv", header=T)
t7 <- read.csv("t7.csv", header=T)
t8 <- read.csv("t8.csv", header=T)
# create one data frame with all rows - all dates and stations
resp2<-rbind(t0, t1, t2, t3, t4, t5, t6, t7, t8)
# Add chamber column 
resp2$chamber<-rep(chamber)
# Add t column
resp2$t <- c(T0, T1, T2, T3, T4, T5, T6, T7,T8)
# Combine both experiments
resp<-rbind(resp1,resp2)
# New column called experiment
ex1 <- rep('ex1_290918', 432)
ex2 <- rep('ex2_011018', 432)
resp$experiment <-c(ex1,ex2)
#make sure chamber, t and experiments are treated as factors
resp$chamber<-as.factor(resp$chamber)
resp$t<-as.factor(resp$t)
resp$experiment <- as.factor(resp$experiment)

resp_mean<-resp%>%
  group_by(experiment,chamber,t)%>%
  summarise(mean_resp=mean(Value))
# add 'real' time column
realtime <- c(0, 3, 6, 9, 12, 15, 18, 21, 24) # hours
resp_mean$time<-rep(realtime, 24)

#Remove values below 100 umol, because we don't want anaerobic processes
resp_mean100<-resp_mean%>%
  filter(mean_resp>=100)

#before plotting: test if regression is significant
resp_mean100<- within(resp_mean100, code3<-(experiment:chamber)[drop = TRUE])
resp.mod1 <-lmList(mean_resp ~ time|code3, na.action = na.omit, data = resp_mean100)
summary(resp.mod1) # All vials show significant decrease in o2 with time


# Get uptake as umol/h by using model coefficients
coef.mod1 <- coef(resp.mod1)
colnames(coef.mod1) <- c("c", "slope") # slope is umol/L/h or uM/h 
# slope is umol/L/h or uM/h and c is the intercept
# we now have a table with slope (m) and b (y-intercept) --> enough to draw a graph

# make slope positive
coef.mod1$slope <- coef.mod1$slope * -1
head(coef.mod1)
#c    slope
#ex1_290918:A01 239.0818 4.376862
#ex1_290918:A02 242.7336 6.353625
#ex1_290918:A03 227.5821 6.301542
#ex1_290918:A04 236.4223 4.539818
#ex1_290918:B05 253.2211 5.678022
#ex1_290918:B06 251.4057 6.368242

#add code in so we know what is what
coef.mod1$code3 <- rownames(coef.mod1)

# remove non-sig slopes --> not necessary as all are significant 
mean(coef.mod1$slope) 
# mean uptake is [1] 6.648551 umol/L/h

# function to calculate standard error of the mean
std <- function(x) sd(x)/sqrt(length(x))

std(coef.mod1$slope) 
#[1] 0.5427785


################# Is there a difference in O2 uptake between experiments? (absolut concentration)
#t.test
Slopes1<- coef.mod1[1:12,2]
Slopes2<- coef.mod1[13:24,2]

t.test(Slopes1, Slopes2, paired = F, var.equal = F)

#Welch Two Sample t-test

#data:  Slopes1 and Slopes2
#t = -2.9246, df = 14.664, p-value = 0.01066
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -4.7661083 -0.7430212
#sample estimates:
#  mean of x mean of y 
#5.271269  8.025833
sd(Slopes1)
#[1] 1.248236
sd(Slopes2)
#[1] 3.014497

# There is a significant difference in oxygen uptake per liter and hour between the batches

# compare O2 decrease to copepod biomass 
cop<- read.csv("Copepods.csv", header=T)
#add in coeff O2 uptakes
cop$slopes <- coef.mod1$slope
#add column of added copepod biomass
cop$weight1 <- 0.46964*cop$Copepod.1-1.12306
cop$weight2<- 0.46964*cop$Copepod.2-1.12306
cop$weighttotal <- cop$weight1+cop$weight2

#install.packages("basicTrendline")
library(basicTrendline)
plot(cop$weighttotal, cop$slopes)
trendline(cop$weighttotal, cop$slopes, model="line2P", show.Rpvalue=T)
#R2 0.04, p=0.34
#no relationship between oxygen decline and biomass (calculated from prosome length)

#compare o2 decrease to c content
cop$carbon1<-41.638*cop$Copepod.1+112.930
cop$carbon2<-41.638*cop$Copepod.2+112.930
cop$carbontotal<- cop$carbon1+cop$carbon2

plot(cop$carbontotal, cop$slopes)
trendline(cop$carbontotal, cop$slopes, model="line2P", show.Rpvalue=T)
# R2 0.04, p=0.34
#no relationship between oxygen decline and carbon content (calculated from prosome lenght)
# both makes sense, as both parameter are derived from the same thing

#compare carbon of different experiments
carbon1<- cop[1:12,14:15]
carbon2<- cop[13:24,14:15]

carbon11<-c(carbon1$carbon1,carbon1$carbon2)
carbon22<-c(carbon2$carbon1,carbon2$carbon2)
t.test(carbon11, carbon22, paired = F, var.equal = F)
#not different!

#compare prosome length of different experiments 
length1<- cop[1:12,c(5,8)]
length2<- cop[13:24,c(5,8)]

length11<-c(length1$Copepod.1, length1$Copepod.2)
length22<-c(length2$Copepod.1, length2$Copepod.2)
t.test(length11, length22, paired=F, var.equal = F)
#not different

#same content in terms of carbon and copepods, but different O2 uptake
#used to plot

resp_mean100$time <- as.factor(resp_mean100$time)
resp_mean100 <- within(resp_mean100, code3<-(experiment:time)[drop = TRUE])

#Take mean of oxygen concentration for each t and measurement time
o2 <- aggregate(resp_mean100$o2, by=list(resp_mean100$code3), mean)
std <- function(x) sd(x)/sqrt(length(x)) # function for std error of the mean
o2.se <- aggregate(resp_mean100$o2, by=list(resp_mean100$code3), std) # std error o2 per temp

# normalise to first o2 concentration for each temperature to compare on plot
o2$norm <- c((o2$x[1:9]/o2$x[1]), (o2$x[10:18]/o2$x[10]))


# 'normalise' std error
o2.se$norm <- o2.se$x/100

#add columns back in
o2$time <- rep(c(0,3,6,9,12,15,18,21,24), 2)
o2$experiment <- c(rep('29/09/18', 9), rep('01/10/18', 9))
o2.se$time <- rep(c(0,3,6,9,12,15,18,21,24), 2)
o2.se$experiment <- c(rep('29/09/18', 9), rep('01/10/18', 9))


#################plotting

Fig2<-plot_grid(plotA_final, plotB_final, plotC_final, labels= 'AUTO')

Fig2_extra
#Renaming column
names(o2)[names(o2)=="experiment"]<-"Experiment"
o2$sd <-o2.se$norm

plotD<-ggplot(data=o2, aes(x=time, y=norm, col=Experiment))+
  geom_point()+
  geom_line()+
  geom_errorbar(aes(ymin=norm-sd, ymax=norm+sd), width=.2, position=position_dodge(0.05))+
  scale_color_manual(values=c("seagreen4","royalblue4"))+
  labs(x='Time (hours)', y=expression(bold(paste('Normalised ',O[2],' concentration'))))+
  theme_cowplot()+
  theme(text=element_text(family="Times New Roman", face="bold", size=14))+
  theme(legend.position=c(0.05,0.2))

Fig2_extra<-plot_grid(plotA_final, plotB_final, plotC_final, plotD, labels= 'AUTO')
Fig2_extra
ggsave("Figure2extra.png", Fig2_extra, width = 1280/72, height = 720/72, dpi = 72)



#calculating the oxygen consumption from umol/L and h to umol/h by multiplying with volume
volume<- c(0.01,0.02,0.02, 0.02, 0.01, 0.02, 0.02, 0.02, 0.01, 0.02, 0.02, 0.02)
coef.mod1$volume <-volume
coef.mod1$consumption<- coef.mod1$slope*coef.mod1$volume
#now we have O2 consumption in umol/h
mean(coef.mod1$consumption)
#[1] 0.1151117 umol h-1
sd(coef.mod1$consumption)
#[1] 0.04440586

mean(coef.mod1$consumption[1:12])
mean(coef.mod1$consumption[13:24])
#[1] 0.095036
#[1] 0.1347879
sd(coef.mod1$consumption[1:12])
sd(coef.mod1$consumption[13:24])

coef.mod1$copug <- cop$carbontotal
coef.mod1$copumol<-coef.mod1$copug/12.0107

coef.mod1$k<-coef.mod1$consumption/coef.mod1$copumol*24

mean(coef.mod1$k)
#[1] 0.06915749
sd(coef.mod1$k)
#[1] 0.02595325
k1<-coef.mod1$k[1:12]
k2<-coef.mod1$k[13:24]

t.test(k1, k2, paired=F, var.equal = F)
#data:  k1 and k2
#t = -2.5006, df = 21.637, p-value = 0.02049
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  -0.043749951 -0.004059664
#sample estimates:
#  mean of x  mean of y 
#0.05720509 0.08110989

#microbial turnover rates k (d-1)
mean(coef.mod1$k[1:12])
#[1] 0.05720509
sd(coef.mod1$k[1:12])
#[1] 0.02184795
max(coef.mod1$k [1:12])
#[1] 0.07914303
min(coef.mod1$k[1:12])
#[1] 0.01466979

mean(coef.mod1$k[13:24])
#[1] 0.08110989
sd(coef.mod1$k[13:24])
#[1] 0.02488671
max(coef.mod1$k [13:24])
#[1] 0.1134663
min(coef.mod1$k[13:24])
#[1] 0.03675046

mean(coef.mod1$k)
#[1] 0.06915749
sd(coef.mod1$k)
#[1] 0.02595325



####Figure 5 - Deep POC flux####
ncin <- nc_open("IMOS_ABOS-SOTS_KF_20180322_SAZ47_FV01_SAZ47-20-2018_PARFLUX-Mark78H-21_END-20180322_C-20200416.nc")
#Get data
Time <- ncvar_get(ncin, "TIME")#Time is number of days since 1950-01-01T00:00:00 UTC
Time_new<-as.POSIXct(as.Date(Time, origin="1950-01-01"))
Depth <- ncvar_get(ncin, "NOMINAL_DEPTH")
POC_flux <- ncvar_get(ncin, "POC_mass_flux")
POC_uncert <- ncvar_get(ncin, "POC_mass_flux_uncertainty")
names(ncin$var)

SOTS <- data.frame(Time_new, Depth, POC_flux, POC_uncert)
SOTS$Depth <-as.factor(Depth)

#Figure 5 plot
SOTS%>%
  filter(Time_new>="2018-06-15"&Time_new<= '2018-10-29')%>%
  ggplot(aes(x=Time_new, y=POC_flux, group=Depth, color=Depth))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=POC_flux-POC_uncert, ymax= POC_flux+POC_uncert), width=.5,
                position=position_dodge(0.05))+
  labs(x="Time", y=bquote('POC flux ('*mg~m^-2~d^-1*')'))+
  scale_colour_manual(values=c("coral", "darkslategray4", "darkslategray"), 
                      labels =c("1000 m","2000 m", "3800 m"))+
  theme_cowplot()+
  theme(legend.position = c(0.05,0.8),axis.title = element_text(size=20),
        axis.text = element_text(size=20), text=element_text(size=20))

#zoom in at the time during the EAC voyage
SOTS%>%
  filter(Time_new>="2018-09-25"&Time_new<= '2018-10-29')%>%
  ggplot(aes(x=Time_new, y=POC_flux, group=Depth, color=Depth))+
  geom_line(size=1)+
  geom_point(size=2)+
  geom_errorbar(aes(ymin=POC_flux-POC_uncert, ymax= POC_flux+POC_uncert), width=.5,
                position=position_dodge(0.05))+
  labs(x="Time", y=bquote('POC flux ('*mg~m^-2~d^-1*')'))+
  scale_colour_manual(values=c("coral", "darkslategray4", "darkslategray"), 
                      labels =c("1000 m","2000 m", "3800 m"))+
  theme_cowplot()+
  theme(legend.position = c(0.05,0.8),axis.title = element_text(size=16),
        text=element_text(size=16))



