check_raster <- function(image, image_dir){
  
  xx <- try(terra::rast(paste0(image_dir,"/",image)))
  
  if(class(xx)[1] == "try-error"){
    return(image)
  } else {
    
    xx <- try(terra::values(terra::rast(paste0(image_dir,"/",image))[[1]]))
    
    if(class(xx)[1] == "try-error"){
      return(image)
    } else {
      return(NULL)
    }
  }
}