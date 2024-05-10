library(lidR)
library(future)
library(EBImage)
library(sf)
library(terra)

borough_name_list <- c("manhattan", "staten", "brooklyn", "queens", "bronx")
#borough_name_list <- c("brooklyn", "queens", "bronx")

for (borough_name in borough_name_list){
  print(borough_name)
  filepath <- paste("/Users/dlm356/dlm356_files/nyc_trees/nyc_lidar_2021/", borough_name, "_test/", sep = "")
  
  setwd(filepath)
  ctg <- readLAScatalog(filepath)
  
  opt_output_files(ctg) <- paste0(tempdir(), "/Users/dlm356/dlm356_files/nyc_trees/nyc_lidar_2021/temp/{ID}") # temp output dir
  opt_chunk_buffer <- 328 # ft == 100 m, not sure if this is actually showing up
  
  las <- readLAS(ctg, select = "xyzrc", filter = "-drop_z_below 0") # this pulls everything into memory...
  
  # planet
  planet_raster <- terra::rast("/Users/dlm356/dlm356_files/nyc_trees/NYC May 2023 test_psscene_analytic_8b_sr_udm2/composite.tif")
  
  # reproject lidar to planet
  las_reproj <- st_transform(las, crs = crs(planet_raster)) # Z is still in US feet (*not* US survey feet)
  las_reproj_z2m <- las_reproj
  las_reproj_z2m$Z <- las_reproj_z2m$Z * 0.3048
  
  pts_classes <- c(LASGROUND, LASWATER, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION)
  las2 <- filter_poi(las_reproj_z2m, Classification %in% c(pts_classes))
  
  # remove unnneeded files and memory
  rm(las, las_reproj)
  gc()
  
  # make DSM
  grid_size <- 0.5 # converting to meters from feet, pretty close to what it was
  dsm <- rasterize_canopy(las2, res = grid_size, pitfree(c(0, 2, 5, 10, 15))) # default, this might take longer than some other methods, overkill?
  
  # make DTM
  dtm <- rasterize_terrain(las2, grid_size, tin())
  
  # Mask out buildings
  mask_base <- dsm*0
  pts_classes <- c(LASBUILDING)
  las_buildings <- filter_poi(las_reproj_z2m, Classification %in% c(pts_classes)) # still needed las_reproj_z2m...
  building_density <- rasterize_density(las_buildings, grid_size) %>% extend(mask_base, fill = 0) %>% crop(mask_base) # extend and then crop to force exact match with the extent of the ground
  mask_building <- mask(mask_base, building_density, inverse = TRUE, maskvalues = 0, updatevalue = 1)
  mask_building_smooth <- focal(x = mask_building, w = 7, fun = "modal") # use 'focal' fn in terra for smoothing, 11 for 0.25, use 7 for 0.5
  mask_building_smooth[is.na(mask_building_smooth)] <- 0 # remove NAs from the mask
  
  dsm_masked <- mask(dsm, mask_building_smooth, maskvalues = 1, updatevalue = 0)
  dtm_masked <- mask(dtm, mask_building_smooth, inverse = TRUE, maskvalues = 1, updatevalue = 0)
  
  #new_dsm <- dsm_masked + dtm_masked
  
  # Inverted density mask for tree CHM to remove everything else
  mask_base <- dsm*0
  pts_classes <- c(LASGROUND, LASWATER, LASHIGHVEGETATION, LASLOWVEGETATION, LASMEDIUMVEGETATION, LASHIGHVEGETATION, LASBUILDING)
  las_fullcls <- filter_poi(las_reproj_z2m, Classification %in% c(pts_classes)) # still needed las_reproj_z2m...
  fullcls_density <- rasterize_density(las_fullcls, grid_size) %>% extend(mask_base, fill = 0) %>% crop(mask_base) # extend and then crop to force exact match with the extent of the ground
  mask_fullcls <- mask(mask_base, fullcls_density, inverse = TRUE, maskvalues = 0, updatevalue = 1)
  mask_fullcls_smooth <- focal(mask_fullcls, w = 7, fun = "modal", na.rm = T) # need to remove NAs otherwise get square holes
  mask_fullcls_smooth[is.na(mask_fullcls_smooth)] <- 0 # remove NAs from the mask
  dsm_masked <- mask(dsm_masked, mask_fullcls_smooth, maskvalues = 0, updatevalue = 0) # invert mask
  dtm_masked2 <- mask(dtm, mask_fullcls_smooth, maskvalues = 1, updatevalue = 0) # this is the fullcls version of dtm, mask again with buildings
  dtm_masked3 <- mask(dtm_masked2, mask_building_smooth, maskvalues = 1, updatevalue = 0)
  dtm_masked_sum <- dtm_masked + dtm_masked3
  new_dsm <- dsm_masked + dtm_masked_sum
  
  # Make CHM
  chm <- new_dsm - dtm # try with new one; 0.25 m might be too small for this point density, might need to do 0.5 m do be safe
  chm[chm < 0] <- 0 # remove negative values
  # note: NA removal causes zero values to be introduced at the borders due to reprojection shift in rect shape of file
  
  # apply smoothing filter on the original CHM because otherwise holes show up
  ker <- matrix(1,5,5)
  chm_s <- focal(chm, w = ker, fun = median, na.rm = T) # need to remove NAs otherwise get square holes
  
  # Do watershed segmentation
  crowns_watershed <- lidR::watershed(chm_s, th_tree = 2, tol = 0.5, ext = 2)() # generates a raster, need to tune this but getting OK results now
  
  # convert to polygon
  poly_crowns_watershed <- terra::as.polygons(crowns_watershed)
  
  # Fill holes
  poly_crowns_watershed_filled <- terra::fillHoles(poly_crowns_watershed)
  
  # write out polygons
  terra::writeVector(poly_crowns_watershed_filled, # write out the FILLED versions
                     filename = paste(borough_name, "_9tiles_test", sep = ""),
                     filetype = "ESRI Shapefile")
  
  # write out CHM, not kernel smoothed
  terra::writeRaster(chm,
                     filename = paste(borough_name, "_9tiles_test.tif", sep = ""))
  
  # clear all variables to free up memory
  rm(list = ls())
  gc()
}