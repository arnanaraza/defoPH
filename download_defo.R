### REMOTE SENSING INPUTS DOWNLOADER USING RGEE

# prelims
rm(list=ls())
pacman::p_load(GADMTools,rgee,sf,raster,rgdal,gdalUtils)
ee_Initialize(email='arnanaraza2006@gmail.com',drive = TRUE)

# globals
main.dir <- 'C:/defoPH/'
setwd(main.dir)

## Funciton to download RADD alerts 
DL <- function(aoi, outdir, rsl){

  # Use block polys Define a region of interest with sf
  aoi_sf <- st_as_sf(aoi)
  st_crs(aoi_sf) <-  4326
  ee_roi <- aoi_sf %>%
    st_geometry() %>%
    sf_as_ee()
  
  # Get alerts
  sat.col <- ee$ImageCollection("projects/radar-wur/raddalert/v1")$
    filterMetadata('layer','contains','alert')$
    filterMetadata('geography','contains','asia')$
    sort("system:time_end", T)$first()
  
  img <-  sat.col$select('Date')
  
  # Block-level download using ee_as_raster, chunks of tiles
  dir.create(file.path(main.dir, outdir))
  outdir <- paste0(main.dir,outdir)
  
  AGB.list <- list()

  
  for (i in 1:length(aoi)){
    dsn <-paste0(outdir,'RADD_', i,'.tif')
    print(dsn)
    a <- aoi[i,]@bbox
    g <- c(a[1],a[2],a[3],a[4])
    geometry <- ee$Geometry$Rectangle(
      coords = g,
      proj = "EPSG:4326",
      geodesic = FALSE)
    
    ee_raster <- ee_as_raster(
      image = img,
      region = geometry,
      dsn = dsn,
      scale = 30,
      via = "drive",
      maxPixels=100000000000)
    
    AGB.list <- c(AGB.list, ee_raster)
  }
  return(AGB.list)
}


#get RADD alerts per AOI
#PH <- readOGR(getwd(), 'aoi')
#PH <- PH[PH$NAME_1 == 'Palawan',]
PH = gadm_sp_loadCountries("PHL", level=1, basefile="./")
poi <- c('Isabela', 'Palawan', 'Quirino' ,'Quezon', 'Laguna', 'Cagayan', 'Apayao',
         'Aurora', 'Nueva Vizcaya', 'Bulacan', 'Rizal')
PH1 <- PH$spdf
aoi <- PH1[PH1$NAME_1 %in% poi, ]
DL(aoi, 'data/RADD/', 30)

dir_radd <- paste0(main.dir, '/data/RADD')
setwd(dir_radd)
r_files <- list.files(getwd(),'tif')
r_files <- r_files[!grepl('vrt|aux|vat', r_files)]
fname <- 'radd_ph.vrt'
gdalbuildvrt(gdalfile = r_files,output.vrt = fname, overwrite=T)


## Funciton to download forest baseline
DL1 <- function(aoi, outdir, rsl){
  
  # Use block polys Define a region of interest with sf
  aoi_sf <- st_as_sf(aoi)
  st_crs(aoi_sf) <-  4326
  ee_roi <- aoi_sf %>%
    st_geometry() %>%
    sf_as_ee()
  
  # Get baseline
  for.base <- ee$ImageCollection("projects/radar-wur/raddalert/v1")$
    filterMetadata('layer','contains','forest_baseline')$
    filterMetadata('geography','contains','asia')$first()
  
  
  # Block-level download using ee_as_raster, chunks of tiles
  dir.create(file.path(main.dir, outdir))
  outdir <- paste0(main.dir,outdir)
  
  AGB.list <- list()
  
  
  for (i in 1:length(aoi)){
    dsn <-paste0(outdir,'FMASK_', i,'.tif')
    a <- aoi[i,]@bbox
    g <- c(a[1],a[2],a[3],a[4])
    geometry <- ee$Geometry$Rectangle(
      coords = g,
      proj = "EPSG:4326",
      geodesic = FALSE)
    
    ee_raster <- ee_as_raster(
      image = for.base,
      region = geometry,
      dsn = dsn,
      scale = rsl,
      via = "drive",
      maxPixels=100000000000)
    
    AGB.list <- c(AGB.list, ee_raster)
  }
  return(AGB.list)
}
DL1(aoi, 'data/ForestMask/', 100)
dir_fmask <- paste0(main.dir, '/data/ForestMask')
setwd(dir_fmask)
f_files <- list.files(getwd(),'tif')
f_files <- f_files[!grepl('vrt|aux|vat', f_files)]
fname <- 'forest_mask_2019.vrt'
gdalbuildvrt(gdalfile = f_files,output.vrt = fname, overwrite=T)
fmask <- raster(fname)
fmask[fmask<=0] <- NA
setwd(paste0(main.dir,'/data'))
writeRaster(fmask, 'forest_mask_2019.tif')

#source...
r <- raster(fname)
r <- mask(r,aoi)
Ras2Pt(r, 'RADD')
