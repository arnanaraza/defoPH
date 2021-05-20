###### WORKFLOW TO PRODUCE FIGURES OF ARAZA ET AL. DEFO HOTSPOT PAPER ########


### 1. Preliminaries

pacman::p_load(raster,rgdal,gdalUtils, plyr,dplyr,ggplot2,lubridate,ggpubr,rgeos)
SRS <- CRS('+init=epsg:4326') 
setwd('C:/defoPH/data')
aoi <- readOGR(getwd(), 'forest_mask1')


### 2. Create forest mask as 30% tree cover and above for hotspot mapping

tc_files <- c('10N_110E_treecover2010_v3.tif', '20N_110E_treecover2010_v3.tif',
              '20N_120E_treecover2010_v3.tif') #pre-downloaded tiles
fname <- 'fmask.vrt'
gdalbuildvrt(gdalfile = tc_files,output.vrt = fname, overwrite=T)
tc <- raster(fname)
aoi <- readOGR(getwd(), 'aoi')
aoi <- aoi[aoi$NAME_1 != 'Ilocos Norte',]
tc_aoi <- crop(tc, extent(aoi))
tc_aoi <- aggregate(tc_aoi, 10, 'mean')
tc_aoi1 <- mask(tc_aoi, aoi)
tc_aoi1[tc_aoi1<30] <- NA
tc_pol <- rasterToPolygons(tc_aoi1, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=T)
tc_pol <- spTransform(tc_pol,  CRS('+init=epsg:32651') )
writeOGR(tc_pol,getwd(), 'forest_mask', 'ESRI Shapefile', overwrite_layer = T)


### 3. GLAD and GFC data VRT (manually downloaded 3 tiles overlapping the study area at http://glad-forest-alert.appspot.com/ and 
### https://earthenginepartners.appspot.com/science-2013-global-forest/download_v1.7.html, respectively)

gfc_files <- list.files(getwd(), 'Hansen')[1:3]
fname <- 'fl_ph.vrt'
gdalbuildvrt(gdalfile = gfc_files,output.vrt = fname, overwrite=T)

glad_files <- list.files(getwd(), 'GLADalert')[7:9]
fname <- 'glad_ph.vrt'
gdalbuildvrt(gdalfile = glad_files,output.vrt = fname, overwrite=T)
glad_files1 <- list.files(getwd(), 'GLADalert')[10:12]
fname1 <- 'glad_ph1.vrt'
gdalbuildvrt(gdalfile = glad_files1,output.vrt = fname1, overwrite=T)


### 4. Create SDPF (points) for the defo alerts (include low confidence alerts (25% of total for RADD))

Ras2Pt <- function(r=r, prod='RADD'){
  
  if  (prod=='RADD'){
    r.pts <- rasterToPoints(r,spatial=T,fun=function(x){x>0})
    r.pts$radd_ph <- r.pts$radd_ph +2000000
    r.pts$yr <- substr(r.pts$radd_ph,1,4)
    r.pts$day <- substr(r.pts$radd_ph,5,7)
    date <- with(r.pts, paste(r.pts$yr, r.pts$day))
    r.pts$date <- strptime(date, "%Y %j")
    r.pts$date1 <- format(r.pts$date, "%m/%d/%Y")
    pf <- readOGR(getwd(), 'forest_mask_2019b')
    r.pts <- raster::intersect(r.pts, pf)
    r.pts=r.pts[,c(1:5)]
    }
  else if (prod=='FL'){
    r.pts <- rasterToPoints(r,spatial=T,fun=function(x){x>0})
    r.pts$date <- as.Date(paste0(r.pts$fl_ph + 2000, '-01-01'))
    r.pts$date1 <- format(r.pts$date, "%m/%d/%Y")
    
  }else{
    r.pts <- rasterToPoints(r[[1]],spatial=T,fun=function(x){x>0})
    names(r.pts) <- 'day'
    r.pts$yr <- '2020'
    date <- with(r.pts, paste(r.pts$yr, r.pts$day))
    r.pts$date <- strptime(date, "%Y %j")
    r.pts$date1 <- format(r.pts$date, "%m/%d/%Y")
    r.pts2020 <- r.pts
    r.pts <- rasterToPoints(r[[2]],spatial=T,fun=function(x){x>0})
    names(r.pts) <- 'day'
    r.pts$yr <- '2021'
    date <- with(r.pts, paste(r.pts$yr, r.pts$day))
    r.pts$date <- strptime(date, "%Y %j")
    r.pts$date1 <- format(r.pts$date, "%m/%d/%Y")
    r.pts <- rbind(r.pts, r.pts2020)
    pf <- readOGR(getwd(), 'forest_mask_2019b')
    r.pts <- raster::intersect(r.pts, pf)
    r.pts=r.pts[,c(1:4)]
  }
  
  #projection to UTM Zone 51N
  r.pts1 <- spTransform(r.pts,  CRS('+init=epsg:32651'))
  aoi <- readOGR(getwd(), 'forest_mask1')
  r.pts2 <- raster::intersect(r.pts1,aoi)
  writeOGR(r.pts2, getwd(), paste0(prod,'_pts'), 'ESRI Shapefile', overwrite_layer = T)
}

Ras2Pt(raster('fl_ph.tif', 'FL'))
Ras2Pt(stack(c('glad_ph.vrt', 'glad_ph1.vrt')) , 'GLAD')
Ras2Pt(raster('./RADD_doy/radd_ph.vrt'), 'RADD')

### 4. Create point data of quarterly deforestation pixels from 2000-2020 separately for RADD and GLAD

Quarterly <- function(prod='GLAD'){
  r.pts=readOGR(getwd(),'FL_pts')
  
  if(prod=='GLAD'){
    r.pts1=readOGR(getwd(), 'GLAD_pts')
  }else{  r.pts1=readOGR(getwd(), 'RADD_pts')}
  
  r.pts$qrtr <- 0
  
  r.pts1$day <- as.numeric(r.pts1$day)
  r.pts1$qrtr <- ifelse(r.pts1$day <= 366, 4,NA)
  r.pts1$qrtr <- ifelse(r.pts1$day <= 273, 3, r.pts1$qrtr )
  r.pts1$qrtr <- ifelse(r.pts1$day <= 181, 2, r.pts1$qrtr )
  r.pts1$qrtr <- ifelse(r.pts1$day <= 90, 1,  r.pts1$qrtr)
  r.pts1$qrtr <- ifelse(r.pts1$yr ==  "2021", 5, r.pts1$qrtr)
  
  fl.alrt <- rbind(r.pts[,c('date1', 'qrtr')], r.pts1[,c('date1', 'qrtr')])
  
  fl.alrt1 <- lapply(sort(unique(fl.alrt$qrtr)), function(x) subset(fl.alrt, fl.alrt$qrtr==x))  
  
  q1 <- rbind(fl.alrt1[[1]], fl.alrt1[[2]])
  print(paste(nrow(q1), 'pts for Q1'))
  
  q2 <- rbind(q1, fl.alrt1[[3]])
  print(paste(nrow(q2), 'pts for Q2'))
  
  q3 <- rbind(q2, fl.alrt1[[4]])
  print(paste(nrow(q3), 'pts for Q3'))
  
  q4 <- rbind(q3, fl.alrt1[[5]])
  print(paste(nrow(q4), 'pts for Q4'))
  
  q5 <- rbind(q4, fl.alrt1[[6]])
  print(paste(nrow(q5), 'pts for Q5'))
  
  fl.alrt2 <- list(q1,q2,q3,q4,q5)
  lapply(1:5, function(x) 
    writeOGR(fl.alrt2[[x]], getwd(), paste0(prod,'_',x), 'ESRI Shapefile',overwrite=T))  
  
}

Quarterly('RADD')
Quarterly('GLAD')



### 5. Create nice cumulative graphs of the alerts (Figure 4 of paper)
setwd('C:/defoPH/data')
pts_list <- c('RADD_pts', 'GLAD_pts')
pts.all <- lapply(1:2, function(x) readOGR(getwd(), pts_list[[x]]))
aoi <- readOGR(getwd(), 'forest_mask1')

a=pts.all[[1]]
a$date <- as.Date(a$date)
a$mo <- month(a$date)
a$defo <- 1
a$id <- paste0(a$yr,'_',a$mo)

a1 <- ddply(as.data.frame(a), .(date), summarize, sum=sum(defo))
a2 <- ddply(as.data.frame(a), .(id), summarize, sum=sum(defo))
a1$sum1 <- cumsum(a1$sum)      
a2$sum1 <- cumsum(a2$sum)      

b=pts.all[[2]]
b <- raster::intersect(b,aoi)
b$date <- as.Date(b$date)
b$mo <- month(b$date)
b$defo <- 1
b$id <- paste0(b$yr,'_',b$mo)

b1 <- ddply(as.data.frame(b), .(date), summarize, sum=sum(defo))
b2 <- ddply(as.data.frame(b), .(id), summarize, sum=sum(defo))
b1$sum1 <- cumsum(b1$sum)   
b2$sum1 <- cumsum(b2$sum)

a1$product <- 'RADD'
b1$product <- 'GLAD'
a2$product <- 'RADD'
b2$product <- 'GLAD'

c1 <- rbind(a1,b1)
c2 <- rbind(a2,b2)

c1 <- c1[with(c1, order(date)), ]
c2 <- c2[with(c2, order(id)), ]

c1$ha <- (c1$sum * 900) / 10000
c1$ha1 <- (c1$sum1 * 900) / 10000

c2$ha <- (c2$sum * 900) / 10000
c2$ha1 <- (c2$sum1 * 900) / 10000

p1 = ggplot(c1,aes(date, ha1,color=product)) +scale_color_manual(values=c("dark green", "dark blue")) +
  geom_line() + theme_bw() +  xlab("Date") + ylab("Cummulative deforestation (ha)") +
  theme(strip.text.x = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),text = element_text(size=16),legend.title = element_blank(),
        legend.position=c(.9,.15)) + scale_x_date(date_breaks = "2 months", date_labels =  "%b %Y") 

c2 <- subset(c2,c2$id != '2020_NA')
p2 = ggplot(c2,aes(id, ha,colour=product))+scale_color_manual(values=c("dark green", "dark blue")) +
  geom_bar(stat='identity', position='dodge',fill='white', width=0.75)+ theme_bw() +  xlab("Year and month") + ylab("Deforested area (ha)") +
  scale_x_discrete(limits=c('2020_1','2020_2','2020_3','2020_4','2020_5','2020_6','2020_7',
                            '2020_8','2020_9','2020_10','2020_11','2020_12','2021_1','2021_2',
                            '2021_3','2021_4')) +  
  theme(strip.text.x = element_blank(),text = element_text(size=16),legend.title = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),legend.position=c(.9,.75))

p <- ggarrange(p1,p2, labels = c("a", "b"),
                    ncol = 1, nrow = 2)
setwd('C:/defoPH/results')
ggsave(p, file=paste0("Figure3a.png"), width = 10, height = 8.5, units = "in")


###################################################################################################################
############### AT THIS POINT THE EMERGING HOTSPOT ANALYSIS MODEL SHOULD BE RUN AT ARCMAP >10.5 ################### 
###################################################################################################################
                  

### 6. Quarterly graph of deforestation pixels inside and outside hotspots (Figure 5 of paper)

a <- subset(a, a$yr == '2020')
b <- subset(b, b$yr == '2020')
a$qrtr <- quarter(a$date)
b$qrtr <- quarter(b$date)

setwd('C:/defoPH/intermediate')
eha.files <- list.files(getwd(), glob2rx("*yr*shp*"))
eha.files0 <- list.files(getwd(), glob2rx("*yr*shp**xml"))
eha.files <- setdiff(eha.files, eha.files0)
eha.files <-  ldply(lapply(eha.files, function(x) tools::file_path_sans_ext(x)),data.frame)
eha.files <- eha.files[-c(5,10),]
eha.all <- lapply(eha.files, function(x) readOGR(getwd(),x))
eha.all1 <- lapply(eha.all, function(x) {x$n <- 1;return(x)})

#defo inside hotspots per quarter
Qrtr <- function(list=eha.all1, pattern='No Pattern Detected'){
  
  eha.qrtr <- lapply(eha.all1, function(x) subset(x, x$PATTERN == pattern))
  eha.defo.glad <-  lapply(1:4, function(x) over(subset(b, b$qrtr == x),eha.qrtr[[x]]))
  eha.defo.radd <-  lapply(5:8, function(x) over(subset(a, a$qrtr == x-4),eha.qrtr[[x]]))
  
  in.defo.glad <-  lapply(eha.defo.glad, function(x) subset(x, !is.na(x$PATTERN)))
  in.defo.radd <-  lapply(eha.defo.radd, function(x) subset(x, !is.na(x$PATTERN)))
  
  gl <- c(sum(in.defo.glad[[1]]$n),sum(in.defo.glad[[2]]$n),sum(in.defo.glad[[3]]$n),
          sum(in.defo.glad[[4]]$n))
  ra <- c(sum(in.defo.radd[[1]]$n),sum(in.defo.radd[[2]]$n),sum(in.defo.radd[[3]]$n),
          sum(in.defo.radd[[4]]$n))
  defo.glad <- data.frame(c(1:4),gl, 'GLAD')
  defo.radd <-data.frame(c(1:4),ra, 'RADD')
  defo.glad$sum1 <- cumsum(defo.glad[[2]])
  defo.radd$sum1 <- cumsum(defo.radd[[2]])
  names(defo.glad) <- names(defo.radd)
  d.inside <- rbind(defo.glad,defo.radd)
  d.inside$ha <-  (d.inside$sum1 * 900) / 10000
  
  names(d.inside) <- c('qrtr' ,'sum','prod', 'sum1', 'ha')
  d.inside$ha1 <-  (d.inside$sum * 900) / 10000
  d.inside
}

patterns <- c("New Hot Spot", "Consecutive Hot Spot", "Sporadic Hot Spot", "No Pattern Detected")
d.all <- lapply(patterns, function(x) Qrtr (eha.all1, x))

p1 = ggplot(data=d.all[[1]],aes(x=qrtr,y=ha1, color=prod)) +scale_color_manual(values=c("dark green", "dark blue")) +
  geom_line(stat='identity',position=position_dodge())+geom_point(size=2)+ theme_bw() + xlab("Quarter, 2020")+
  ylab("Deforested area (ha)") +   ylim(0, 3800)+ ggtitle("New")+
  theme(text = element_text(size=18), panel.grid.major = element_blank(),
        plot.title = element_text(size=18,vjust = - 10,hjust = +.9),
        panel.grid.minor = element_blank(), legend.position='none')
p2 = ggplot(data=d.all[[2]],aes(x=qrtr,y=ha1, color=prod)) +scale_color_manual(values=c("dark green", "dark blue")) +
  geom_line(stat='identity',position=position_dodge())+geom_point(size=2)+ theme_bw() + xlab("Quarter, 2020")+
  ylab("Deforested area (ha)") +   ylim(0, 3800)+ ggtitle("Consecutive")+
  theme(text = element_text(size=18),axis.text.y = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(size=18,vjust = - 10,hjust = +.9),legend.title = element_blank(),
        panel.grid.minor = element_blank(),legend.position=c(.85,.75))
p3 = ggplot(data=d.all[[3]],aes(x=qrtr,y=ha1, color=prod)) +scale_color_manual(values=c("dark green", "dark blue")) +
  geom_line(stat='identity',position=position_dodge())+geom_point(size=2)+ theme_bw() + xlab("Quarter, 2020")+
  ylab("Deforested area (ha)") +   ylim(0, 3800)+ ggtitle("Sporadic")+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), 
        plot.title = element_text(size=18,vjust = - 10,hjust = +.9),
        panel.grid.minor = element_blank(), legend.position='none')
p4 = ggplot(data=d.all[[4]],aes(x=qrtr,y=ha1, color=prod)) +scale_color_manual(values=c("dark green", "dark blue")) +
  geom_line(stat='identity',position=position_dodge())+geom_point(size=2)+ theme_bw() + xlab("Quarter, 2020")+
  ylab("Deforested area (ha)") +   ylim(0, 3800)+ ggtitle("Outside")+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), axis.text.y = element_blank(),
        plot.title = element_text(size=18,vjust = - 10,hjust = +.9),
        panel.grid.minor = element_blank(), legend.position='none')
setwd('C:/defoPH/results')
p <- ggarrange(p1,p2,p3,p4, labels = c("a", "b", 'c', 'd'),ncol = 2, nrow = 2, font.label = list(size = 20))
ggsave(p, file=paste0("Figure4_fin.png"), width = 12, height = 10, units = "in")

d.all[[1]]$PATTERN <- 'New'
d.all[[2]]$PATTERN <- 'Consecutive'
d.all[[3]]$PATTERN <- 'Sporadic'
d.all[[4]]$PATTERN <- 'Outside'
write.csv(rbind(d.all[[1]],d.all[[2]],d.all[[3]],d.all[[4]]), 'eha.all3.csv',row.names=F)
    

### 7. Tally per hotspot class per alert product (Figure 3 and Table 2 of paper)

setwd('C:/defoPH/intermediate')
glad.pt <- readOGR(getwd(), 'eha_glad_fl_4yr')
glad.pt$pt <- gCentroid(glad.pt, byid=TRUE)

#forests 30% TCD
setwd('C:/defoPH/data')
aoi <- readOGR(getwd(), 'forest_mask1')
glad.pt$over.gl <- over(glad.pt,aoi)[[1]]
glad.pt$n <- 1
glad.pt$Id <- ifelse(glad.pt$Id == 0, 1)
glad.pt$Id <- ifelse(is.na(glad.pt$Id),1 ,0)
glad.fmask <- aggregate(n ~ PATTERN, FUN = sum, data = glad.pt, na.rm = TRUE)

#primary 2020
aoi <- aggregate(raster('forest_mask_2019.tif'), 10, fun='mean')
aoi.pri <- rasterToPolygons(aoi, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=T)
aoi.pri <- spTransform(aoi.pri,  CRS('+init=epsg:32651') )
glad.pt$over.pri <- over(glad.pt,aoi.pri)[[1]]
glad.pt$over.pri <- ifelse(is.na(glad.pt$over.pri),0,glad.pt$over.pri)
glad.pri <- aggregate(n ~ PATTERN+over.pri, FUN = sum, data = glad.pt, na.rm = TRUE)

#pa
aoi <- readOGR(getwd(), 'pa_aoi')
glad.pt$over.pa <- over(glad.pt,aoi)
glad.pt$over.pa <- ifelse(is.na(glad.pt$over.pa$PANAME),0,1)
glad.pa <- aggregate(n ~ PATTERN+over.pa, FUN = sum, data = glad.pt, na.rm = TRUE)
glad.pa


#forests 30% TCD
setwd('C:/defoPH/intermediate')
radd.pt <- readOGR(getwd(), 'eha_rad_fl_4yr')
radd.pt$pt <- gCentroid(radd.pt, byid=TRUE)
radd.pt$n <- 1
radd.fmask <- aggregate(n ~ PATTERN, FUN = sum, data = radd.pt, na.rm = TRUE)

#primary 2020
setwd('C:/defoPH/data')
radd.pt$over.pri <- over(radd.pt,aoi.pri)[[1]]
radd.pt$over.pri <- ifelse(is.na(radd.pt$over.pri),0,radd.pt$over.pri)
radd.pri <- aggregate(n ~ PATTERN+over.pri, FUN = sum, data = radd.pt, na.rm = TRUE)

#pa
aoi.pa <- readOGR(getwd(), 'pa_aoi')
radd.pt$over.pa <- over(radd.pt,aoi.pa)
radd.pt$over.pa <- ifelse(is.na(radd.pt$over.pa$PANAME),0,1)
radd.pa <- aggregate(n ~ PATTERN+over.pa, FUN = sum, data = radd.pt, na.rm = TRUE)
radd.pa

## protected areas graph
setwd('C:/defoPH/data')
muni <- readOGR(getwd(), 'aoi_muni')
radd.pt$muni <- over(radd.pt,muni)[[5]]
radd.muni.pa <- aggregate(n ~muni+over.pa+PATTERN , FUN = sum, 
                          data =subset(radd.pt, radd.pt$PATTERN != 'No Pattern Detected'), na.rm = TRUE)
radd.muni.pa$prod <-'RADD'
radd.pt$muni <- over(radd.pt,muni)[[5]]

glad.pt$muni <- over(glad.pt,muni)[[5]]
glad.muni.pa <- aggregate(n ~muni+over.pa+PATTERN , FUN = sum, 
                          data =subset(glad.pt, glad.pt$PATTERN != 'No Pattern Detected'), na.rm = TRUE)
glad.muni.pa$prod <-'GLAD'

over.muni.pa <- rbind(glad.muni.pa,radd.muni.pa)
over.muni.pa$over.pa <-ifelse(over.muni.pa$over.pa == 1, 'PA', 'Non-PA')

p4 = ggplot(over.muni.pa,aes(reorder(muni, n), n,fill=over.pa))+scale_fill_manual(values=c("grey60", "grey90")) +
  geom_bar(stat='identity',position='dodge', width=0.75)+ theme_bw() +  ylab("Hotspot (n)") +xlab("Province")+
  theme(text = element_text(size=18), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        legend.title = element_blank(), legend.position=c(.85,.15))+  coord_flip()+
  facet_wrap(vars(prod),scales = "free_x")

p4
setwd('C:/defoPH/results')
ggsave(p4, file=paste0("Fig3_fin.png"), width = 8, height = 4, units = "in")

