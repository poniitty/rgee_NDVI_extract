# Mask bad pixels, removes images with low cover of good data
# image <- "aland_C2_T1_LT05_190018_19920608_090529GMT.tif"
mask_NDVI <- function(image, image_dir, ndvi_dir, mincover = 2){
  
  require(terra)
  mask_clear <- function(x) {as.numeric(intToBits(x)[7])}
  
  rs <- rast(paste0(image_dir,"/",image))
  rsn <- names(rs)
  
  cmask <- rs[[grepl("qa_",rsn, ignore.case = T)]]
  cmask[] <- unlist(lapply(as.numeric(values(rs[[grepl("qa_",rsn, ignore.case = T)]])), mask_clear))
  
  rs <- ifel(cmask == 0, NA, rs)
  
  clprop <- mean(values(!is.na(rs[[1]])))*100
  
  if(clprop >= mincover){
    rs <- rs[[1]]
    names(rs) <- "NDVI"
    rs[rs <= 0] <- NA
    writeRaster(rs, paste0(ndvi_dir,"/",image), datatype = "INT2U", overwrite = T)
    return("preserved")
  } else {
    unlink(paste0(image_dir,"/",image), force = T)
    return("removed")
  }
}

# i <- "2022-08-29"
#
extract_landsat_NDVI_gee <- function(i, aoi_ee, epsg, drive_folder, area_dl_dir, name){
  # LANDSAT TM4
  if(year(i) %in% 1984:1993){
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT04/C02/T1_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
    e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
    
    if(class(e)[1] != "try-error"){
      maximg <- dataset$reduce(ee$Reducer$max())
      ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
      ndvi <- ndvi$
        multiply(1000)$ # Multiply before transforming to integers
        toUint16()
      maximg <- maximg$addBands(ndvi)$
        select("nd","QA_PIXEL_max")$
        reproject(crs = paste0("EPSG:",epsg), scale = 30)$
        toUint16()
      
      task_img <- ee_image_to_drive(
        image = maximg,
        fileFormat = "GEO_TIFF",folder = drive_folder,
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
      )
      
      task_img$start()
      toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
      if(class(toe)[1] == "try-error"){
        task_img$cancel()
        return("TIMEOUT")
      } else {
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                 format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
        
        return("TM4 found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("TM4 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("TM4 found")
          }
        }
      }
      
    } else {
      # T2, COLLECTION 2
      dataset <- ee$ImageCollection('LANDSAT/LT04/C02/T2_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
      dataset <- dataset$filterBounds(geometry = aoi_ee)
      dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
      dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
      e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
      
      if(class(e)[1] != "try-error"){
        maximg <- dataset$reduce(ee$Reducer$max())
        ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
        ndvi <- ndvi$
          multiply(1000)$ # Multiply before transforming to integers
          toUint16()
        maximg <- maximg$addBands(ndvi)$
          select("nd","QA_PIXEL_max")$
          reproject(crs = paste0("EPSG:",epsg), scale = 30)$
          toUint16()
        
        task_img <- ee_image_to_drive(
          image = maximg,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
        )
        
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("TM4 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("TM4 found")
            if(class(de)[1] == "try-error"){
              de <- try(ee_drive_to_local(task = task_img, 
                                          dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                       tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                       format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
              
              return("TM4 found")
            }
          }
        }
        
      } else {
        
      }
    }
  }
  
  # LANDSAT TM5
  # All same tasks are done for TM5
  if(year(i) %in% 1984:2013){
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LT05/C02/T1_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
    e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
    
    if(class(e)[1] != "try-error"){
      maximg <- dataset$reduce(ee$Reducer$max())
      ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
      ndvi <- ndvi$
        multiply(1000)$ # Multiply before transforming to integers
        toUint16()
      maximg <- maximg$addBands(ndvi)$
        select("nd","QA_PIXEL_max")$
        reproject(crs = paste0("EPSG:",epsg), scale = 30)$
        toUint16()
      
      task_img <- ee_image_to_drive(
        image = maximg,
        fileFormat = "GEO_TIFF",folder = drive_folder,
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
      )
      
      task_img$start()
      toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
      if(class(toe)[1] == "try-error"){
        task_img$cancel()
        return("TIMEOUT")
      } else {
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                 format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
        
        return("TM5 found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("TM5 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("TM5 found")
          }
        }
      }
      
    } else {
      # T2, COLLECTION 2
      dataset <- ee$ImageCollection('LANDSAT/LT05/C02/T2_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
      dataset <- dataset$filterBounds(geometry = aoi_ee)
      dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
      dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
      e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
      
      if(class(e)[1] != "try-error"){
        maximg <- dataset$reduce(ee$Reducer$max())
        ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
        ndvi <- ndvi$
          multiply(1000)$ # Multiply before transforming to integers
          toUint16()
        maximg <- maximg$addBands(ndvi)$
          select("nd","QA_PIXEL_max")$
          reproject(crs = paste0("EPSG:",epsg), scale = 30)$
          toUint16()
        
        task_img <- ee_image_to_drive(
          image = maximg,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
        )
        
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("TM5 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("TM5 found")
            if(class(de)[1] == "try-error"){
              de <- try(ee_drive_to_local(task = task_img, 
                                          dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                       tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                       format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
              
              return("TM5 found")
            }
          }
        }
        
      } else {
        
      }
    }
  }
  
  # LANDSAT ETM7
  if(year(i) %in% 1999:year(Sys.time())){
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LE07/C02/T1_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
    e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
    
    if(class(e)[1] != "try-error"){
      maximg <- dataset$reduce(ee$Reducer$max())
      ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
      ndvi <- ndvi$
        multiply(1000)$ # Multiply before transforming to integers
        toUint16()
      maximg <- maximg$addBands(ndvi)$
        select("nd","QA_PIXEL_max")$
        reproject(crs = paste0("EPSG:",epsg), scale = 30)$
        toUint16()
      
      task_img <- ee_image_to_drive(
        image = maximg,
        fileFormat = "GEO_TIFF",folder = drive_folder,
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
      )
      
      task_img$start()
      toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
      if(class(toe)[1] == "try-error"){
        task_img$cancel()
        return("TIMEOUT")
      } else {
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                 format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
        
        return("ETM7 found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("ETM7 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("ETM7 found")
          }
        }
      }
      
    } else {
      # T2, COLLECTION 2
      dataset <- ee$ImageCollection('LANDSAT/LE07/C02/T2_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
      dataset <- dataset$filterBounds(geometry = aoi_ee)
      dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B7","ST_B6","QA_PIXEL","QA_RADSAT")
      dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
      e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
      
      if(class(e)[1] != "try-error"){
        maximg <- dataset$reduce(ee$Reducer$max())
        ndvi <- maximg$normalizedDifference(c("SR_B4_max","SR_B3_max"))
        ndvi <- ndvi$
          multiply(1000)$ # Multiply before transforming to integers
          toUint16()
        maximg <- maximg$addBands(ndvi)$
          select("nd","QA_PIXEL_max")$
          reproject(crs = paste0("EPSG:",epsg), scale = 30)$
          toUint16()
        
        task_img <- ee_image_to_drive(
          image = maximg,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
        )
        
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("ETM7 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("ETM7 found")
            if(class(de)[1] == "try-error"){
              de <- try(ee_drive_to_local(task = task_img, 
                                          dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                       tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                       format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
              
              return("ETM7 found")
            }
          }
        }
        
      } else {
        
      }
    }
  }
  
  # LANDSAT OLI8
  if(year(i) %in% 2013:year(Sys.time())){
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC08/C02/T1_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
    e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
    
    if(class(e)[1] != "try-error"){
      maximg <- dataset$reduce(ee$Reducer$max())
      ndvi <- maximg$normalizedDifference(c("SR_B5_max","SR_B4_max"))
      ndvi <- ndvi$
        multiply(1000)$ # Multiply before transforming to integers
        toUint16()
      maximg <- maximg$addBands(ndvi)$
        select("nd","QA_PIXEL_max")$
        reproject(crs = paste0("EPSG:",epsg), scale = 30)$
        toUint16()
      
      task_img <- ee_image_to_drive(
        image = maximg,
        fileFormat = "GEO_TIFF", folder = drive_folder,
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
      )
      
      task_img$start()
      toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
      if(class(toe)[1] == "try-error"){
        task_img$cancel()
        return("TIMEOUT")
      } else {
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                 format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
        
        return("OLI8 found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("OLI8 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("OLI8 found")
          }
        }
      }
      
    } else {
      # T2, COLLECTION 2
      dataset <- ee$ImageCollection('LANDSAT/LC08/C02/T2_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
      dataset <- dataset$filterBounds(geometry = aoi_ee)
      dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
      dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
      e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
      
      if(class(e)[1] != "try-error"){
        maximg <- dataset$reduce(ee$Reducer$max())
        ndvi <- maximg$normalizedDifference(c("SR_B5_max","SR_B4_max"))
        ndvi <- ndvi$
          multiply(1000)$ # Multiply before transforming to integers
          toUint16()
        maximg <- maximg$addBands(ndvi)$
          select("nd","QA_PIXEL_max")$
          reproject(crs = paste0("EPSG:",epsg), scale = 30)$
          toUint16()
        
        task_img <- ee_image_to_drive(
          image = maximg,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
        )
        
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("OLI8 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("OLI8 found")
            if(class(de)[1] == "try-error"){
              de <- try(ee_drive_to_local(task = task_img, 
                                          dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                       tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                       format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
              
              return("OLI8 found")
            }
          }
        }
        
      } else {
        return("OLI8 not found")
      }
    }
  }
  
  # LANDSAT OLI9
  if(year(i) %in% 2021:year(Sys.time())){
    # T1, COLLECTION 2
    dataset <- ee$ImageCollection('LANDSAT/LC09/C02/T1_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
    dataset <- dataset$filterBounds(geometry = aoi_ee)
    dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
    dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
    e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
    
    if(class(e)[1] != "try-error"){
      maximg <- dataset$reduce(ee$Reducer$max())
      ndvi <- maximg$normalizedDifference(c("SR_B5_max","SR_B4_max"))
      ndvi <- ndvi$
        multiply(1000)$ # Multiply before transforming to integers
        toUint16()
      maximg <- maximg$addBands(ndvi)$
        select("nd","QA_PIXEL_max")$
        reproject(crs = paste0("EPSG:",epsg), scale = 30)$
        toUint16()
      
      task_img <- ee_image_to_drive(
        image = maximg,
        fileFormat = "GEO_TIFF", folder = drive_folder,
        region = aoi_ee,
        scale = 30,
        fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
      )
      
      task_img$start()
      toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
      if(class(toe)[1] == "try-error"){
        task_img$cancel()
        return("TIMEOUT")
      } else {
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                 format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
        
        return("OLI9 found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("OLI9 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T1_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("OLI9 found")
          }
        }
      }
      
    } else {
      # T2, COLLECTION 2
      dataset <- ee$ImageCollection('LANDSAT/LC09/C02/T2_L2')$filterDate(i, as.character(as.Date(i)+days(1)))
      dataset <- dataset$filterBounds(geometry = aoi_ee)
      dataset <- dataset$select("SR_B1","SR_B2","SR_B3","SR_B4","SR_B5","SR_B6","SR_B7","ST_B10","QA_PIXEL","QA_RADSAT")
      dataset <- dataset$filterMetadata("CLOUD_COVER", "less_than", 80)
      e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
      
      if(class(e)[1] != "try-error"){
        maximg <- dataset$reduce(ee$Reducer$max())
        ndvi <- maximg$normalizedDifference(c("SR_B5_max","SR_B4_max"))
        ndvi <- ndvi$
          multiply(1000)$ # Multiply before transforming to integers
          toUint16()
        maximg <- maximg$addBands(ndvi)$
          select("nd","QA_PIXEL_max")$
          reproject(crs = paste0("EPSG:",epsg), scale = 30)$
          toUint16()
        
        task_img <- ee_image_to_drive(
          image = maximg,
          fileFormat = "GEO_TIFF",folder = drive_folder,
          region = aoi_ee,
          scale = 30,
          fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
        )
        
        task_img$start()
        toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
        if(class(toe)[1] == "try-error"){
          task_img$cancel()
          return("TIMEOUT")
        } else {
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                   format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
          
          return("OLI9 found")
          if(class(de)[1] == "try-error"){
            de <- try(ee_drive_to_local(task = task_img, 
                                        dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                     tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                     format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
            
            return("OLI9 found")
            if(class(de)[1] == "try-error"){
              de <- try(ee_drive_to_local(task = task_img, 
                                          dsn = paste0(area_dl_dir, "/",name,"_C2_T2_",
                                                       tail(str_split(ei$img_id,"/")[[1]],1),"_",
                                                       format(ei$img_time_start, "%H%M%S%Z"))), silent = T)
              
              return("OLI9 found")
            }
          }
        }
        
      } else {
        return("OLI9 not found")
      }
    }
  }
}