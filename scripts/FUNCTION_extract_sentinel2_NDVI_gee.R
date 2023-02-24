# i <- "2021-08-15"
#

S2_clean <- function(img) {
  # Calculate the NDVI
  ndvi_values <- img$normalizedDifference(c("B8","B4"))
  
  # Extract the quality band
  img_qa <- img$select("SCL")
  
  # Create a mask considering: snow , water, cloud shandows, medium&high clouds, cirrus etc
  # See the mask values in: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_S2_SR_HARMONIZED
  cloud_mask <- img_qa$eq(list(1:3, 6:11))$reduce(ee$Reducer$sum())$gt(0)
  
  # Mask pixels with value zero.
  ndvi_values %>%
    ee$Image$updateMask(cloud_mask) %>%
    ee$Image$copyProperties(img, list("system:time_start"))
}

extract_sentinel2_NDVI_gee <- function(i, aoi_ee, epsg, drive_folder, area_dl_dir, name){
  # Sentinel-2
  
  # Surface reflectance, Harmonized
  dataset <- ee$ImageCollection('COPERNICUS/S2_SR_HARMONIZED')$filterDate(i, as.character(as.Date(i)+days(1)))
  dataset <- dataset$filterBounds(geometry = aoi_ee)
  dataset <- dataset$select("B4","B8","SCL")
  dataset <- dataset$filterMetadata("HIGH_PROBA_CLOUDS_PERCENTAGE", "less_than", 80)
  e <- try(withTimeout({ei <- ee_print(dataset, quiet = T)}, timeout = info_timeout), silent = T)
  
  if(class(e)[1] != "try-error"){
    
    maximg <- dataset$
      map(S2_clean)$
      reduce(ee$Reducer$max())$
      reproject(crs = paste0("EPSG:",epsg), scale = 10)$
      multiply(1000)$ # Multiply before transforming to integers
      toUint16()
    
    task_img <- ee_image_to_drive(
      image = maximg,
      fileFormat = "GEO_TIFF",folder = drive_folder,
      region = aoi_ee,
      scale = 10,
      fileNamePrefix = paste0(name,"_",tail(str_split(ei$img_id,"/")[[1]],1))
    )
    
    task_img$start()
    toe <- try(withTimeout(ee_monitoring(task_img, max_attempts = (my_timeout/5)+1), timeout = my_timeout), silent = T)
    if(class(toe)[1] == "try-error"){
      task_img$cancel()
      return("TIMEOUT")
    } else {
      de <- try(ee_drive_to_local(task = task_img, 
                                  dsn = paste0(area_dl_dir, "/",name,"_S2_SR_harm_",
                                               tail(str_split(ei$img_id,"/")[[1]],1))), silent = T)
      
      return("S2 SR scene found")
      if(class(de)[1] == "try-error"){
        de <- try(ee_drive_to_local(task = task_img, 
                                    dsn = paste0(area_dl_dir, "/",name,"_S2_SR_harm_",
                                                 tail(str_split(ei$img_id,"/")[[1]],1))), silent = T)
        
        return("S2 SR scene found")
        if(class(de)[1] == "try-error"){
          de <- try(ee_drive_to_local(task = task_img, 
                                      dsn = paste0(area_dl_dir, "/",name,"_S2_SR_harm_",
                                                   tail(str_split(ei$img_id,"/")[[1]],1))), silent = T)
          
          return("S2 SR scene found")
        }
      }
    }
  } else {
    return("S2 SR data not found on specified date")
  }
}