# Goal: extract all tree data from images while limiting the amount of i/o for image loading

# Load tree polygons
# Run through Planet extents to get full list of images that overlap with the full set of tree polygons
# Loop through by image instead of by polygon
# For each image
# Load image
# Get all tree crowns that overlap within it
# extract pixel values and weights for all tree crowns within image and metadata image
# write this dataframe out to a file with tree ids and image id
# move to next image
# afterwards
# can load as many (or all) of the outputted dataframes as desired
# can merge into bigger dataframe
# can get date from the image id
# can filter by meta information like clear days
# can calculate ndvi and whatever else we want.

#---

# Load libraries
library(tidyverse)
library(terra) # try doing this all in terra
library(future)
library(future.apply)


#---

# Load lidar polygons
# would need to do tree attribution for polygons too
tree_polys <- vect("/Users/dlm356/dlm356_files/nyc_trees/nyc_lidar_2021/bronx_test/bronx_12255_test_dalponte_fitree/bronx_12255_test_dalponte_fitree.shp")

tree_poly_path <- "/Users/dlm356/dlm356_files/nyc_trees/nyc_lidar_2021/bronx_test/bronx_12255_test_dalponte_fitree/bronx_12255_test_dalponte_fitree.shp"

output_path <- "/Volumes/NYC_geo/Planet/tests/test_extract_output7/" # output locations for all files

# Planet data
setwd("/Volumes/NYC_geo/Planet/raw_images/time_series")
#setwd("2019") # /output_20190701_4b/files")

# Full extent of NYC area loaded for planet
nyc_vect <- vect("/Volumes/NYC_geo/Planet/outline/nyc_buffer_dissolved_filled_convex_hull_latlon.geojson") # full extent of NYC area used to get Planet data

# load json for extents
json_files <- list.files(pattern = "\\metadata.json$", recursive = TRUE)

# load all these Planet extent spatvectors to a single file
all_img_ext_poly <- vect(json_files[1]) # this is just the first one.
for (i in 2:length(json_files)){
  all_img_ext_poly <- rbind(all_img_ext_poly, vect(json_files[i])) # loop to add all extents into single file, may be a better way to do this but it works
}

# need to convert crs
# convert the tree polygon because that will be faster than converting every Planet raster
img_ext_poly <- vect(json_files[1])
tree_polys_reproj <- project(tree_polys, crs(img_ext_poly))

tree_polys_ch <- convHull(tree_polys_reproj) # calculate convex hull of tree polys, might replace this with nyc_vect in the end since these will be close to the same thing...

# Filter Planet images to which overlap with the tree crowns (using convex hull)
tree_intersection <- intersect(tree_polys_ch, all_img_ext_poly) # no intersection means they don't overlap, this reduces the list
tree_intersection_df <- as.data.frame(tree_intersection) # this is the df list of images left over that are available
tree_img_id_list <- tree_intersection_df$id # img list, these are the planet images to load

###
# Load list of all Planet image(s) and the json files
# Setup a different labeling scheme for the 4b vs the 8b scenes and keep them separate and specify which is which
tif_files <- list.files(pattern = "\\.tif$", recursive = TRUE)
json_files <- list.files(pattern = "\\metadata.json$", recursive = TRUE) # already loaded this but maybe better to do it this way instead

# these are different lengths and they should not be, should be json_files x2  = tif_files
# json files is slightly longer, but that's OK because we're doing string detect on tifs

# Cannot pass spatvector and spatraster to a future worker, so need to load it within the function
# https://cran.r-project.org/web/packages/future/vignettes/future-4-non-exportable-objects.html

extractionSubRoutine <- function(nb, tif_files_vec, json_files_vec, output_path, file_prefix, i, tree_img_id_list){
  img_file <- tif_files_vec[str_detect(tif_files_vec, "harmonized_clip.tif")] # this pulls both the 4b and 8b versions
  meta_file <- tif_files_vec[str_detect(tif_files_vec, "udm2_clip.tif")] # this will combine and duplicate
  img_with_meta <- c(rast(img_file), rast(meta_file)) # stack the img and the img_meta so only extract once
  img_convhull <- vect(json_files_vec) # only 1 file!
  tree_polys_reproj2 <- project(vect(tree_poly_path), crs(img_with_meta)) # this needs to be loaded within the function
  
  # Planet json extents are different projection than the actual image files
  img_convhull_reproj2 <- project(img_convhull, crs(img_with_meta)) # need to do it for the convex hull extent every time
  
  # get subset for trees that are within the convex hull for the image
  trees_in_img <- intersect(tree_polys_reproj2, img_convhull_reproj2) # no intersection means they don't overlap, this reduces the list
  # what about edges?? and clipped trees?? could be same date within overlapping scenes...
  # after the fact, take the image that has more pixels for a given tree?
  
  # extract reflectance from img and metadata from img_meta for all tree polygons within the extent (convex hull boundary file for the image)
  tree_vals <- extract(img_with_meta, trees_in_img, method = "simple", weights = TRUE)
  tree_vals <- tree_vals[complete.cases(tree_vals),] # remove rows that are NA in reflectance
  # will need to double check that the tree IDs that are getting exported are unique across ALL images (and don't get relabeled per each image)
  #tree_vals$img_id <- tree_img_id_list[i] # might not need to do this if we include the image id in the file name, save a little space (a lot of space, 36%)
  
  # write output
  output_path_filename <- paste(output_path, file_prefix, tree_img_id_list[i], "_", nb, ".csv", sep = "")
  write_csv(tree_vals, output_path_filename) # paste together output path, img_id into filename
}

extractTreesFromPlanetFuture <- function(i, tree_img_id_list, tif_files, json_files, tree_poly_path, output_path){
  # Load the tif files for the img and meta file
  tif_files_sub <- tif_files[str_detect(tif_files, tree_img_id_list[i])]
  json_files_sub <- json_files[str_detect(json_files, tree_img_id_list[i])]
  
  # Do a 4b vs. 8b detection in here
  tif_files_sub_4b <- tif_files_sub[which(str_detect(tif_files_sub, "_4b/"))]
  tif_files_sub_8b <- tif_files_sub[which(str_detect(tif_files_sub, "_8b/"))]
  json_files_sub_4b <- json_files_sub[which(str_detect(json_files_sub, "_4b/"))]
  json_files_sub_8b <- json_files_sub[which(str_detect(json_files_sub, "_8b/"))]
  
  # 4b routine, original recipe
  if (length(tif_files_sub_4b) > 0 & length(json_files_sub_4b) > 0){
    extractionSubRoutine("4b", tif_files_sub_4b, json_files_sub_4b, output_path, "bronx_12255_test_dalponte_fitree_", i, tree_img_id_list)
  }
  
  # 8b routine, new and improved
  if (length(tif_files_sub_8b) > 0 & length(json_files_sub_8b) > 0){
    extractionSubRoutine("8b", tif_files_sub_8b, json_files_sub_8b, output_path, "bronx_12255_test_dalponte_fitree_", i, tree_img_id_list)
  }
}

# Set up the parallel plan
plan(multisession)  # Use multisession plan for parallel execution

# Execute the processing in parallel
#options(future.globals.onReference = "error") # turn this on for debugging
#options(future.globals.onReference = "ignore") 
future_lapply(1:length(tree_img_id_list), function(i) extractTreesFromPlanetFuture(i, tree_img_id_list, tif_files, json_files, tree_poly_path, output_path))

#future_lapply(1:20, function(i) extractTreesFromPlanetFuture(i, tree_img_id_list, tif_files, json_files, tree_poly_path, output_path))
# tomorrow test run for all time for planet for these polygons

#####
# loading and plotting for NDVI from planet

setwd(output_path)
csv_file_list <- list.files(pattern = "_4b\\.csv") #list.files("\\.csv$")

prefix <- "bronx_12255_test_dalponte_fitree"
date_start <- nchar(prefix)+2
date_length <- 8

for(i in 1:length(csv_file_list)){
#for(i in 1:318){  # test for everything before we get to the 8 band imagery too
#for(i in c(1, 318:320)){
  print(i)
  tree_df <- read_csv(csv_file_list[i], show_col_types = FALSE) # show_col_types = FALSE
  if (nrow(tree_df) > 0){ # do this to make sure there is data in the file
    tree_df$output_file <- csv_file_list[i]
    tree_df$date <- ymd(substr(csv_file_list[i], date_start, date_start + date_length - 1))
    if (i == 1){
      all_tree_df <- tree_df
    } else {
      all_tree_df <- bind_rows(all_tree_df, tree_df)
    }
  }
}
# need to make sure we can use the ID in this data frame to tie back to the original tree polygons...

# Need to setup aggregation for weighted mean again, group_by the tree ID and by the image extracted

# calculate NDVI and make plots

# see if we can see which are likely trees and which aren't if we wanted to do a post-hoc greenness knockout routine


# retain only clear flags
all_tree_df_clear <- subset(all_tree_df, subset = all_tree_df$clear == 1)

# get weighted average
# will need to change this to 8 band versions
all_tree_df_clear_sub <- subset(all_tree_df_clear, select = c("ID", "blue", "green", "red", "nir", "date", "weight"))
all_tree_df_clear_sub_wmean <- all_tree_df_clear_sub %>% group_by(ID, date) %>% summarize(across(c("blue", "green", "red", "nir"), ~ weighted.mean(.x, w = weight)), .groups = "keep") # this did it!

# calculate NDVI
getNDVI <- function(nir, red){
  return((nir - red)/(nir + red))
}
all_tree_df_clear_sub_wmean$ndvi <- getNDVI(all_tree_df_clear_sub_wmean$nir, all_tree_df_clear_sub_wmean$red)

df_plot <- subset(all_tree_df_clear_sub_wmean, subset = ID == 3)

new_theme <- theme_set(theme_bw())
theme_update(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# plot
ggplot(df_plot)+
  geom_line(aes(x = date, y = ndvi), col = "forestgreen") +
  geom_point(aes(x = date, y = ndvi), col = "forestgreen") +
  # geom_line(aes(x = date, y = ndvi, color = as.factor(ID))) +
  # geom_point(aes(x = date, y = ndvi, color = as.factor(ID))) +
  labs(x = "Date", y = "NDVI") + 
  scale_x_date(date_breaks = "3 months", date_labels = "%m-%Y", date_minor_breaks = "1 month",
               limits = as.Date(c("2018-01-01","2024-01-01")), expand = expansion(0))
ggsave("../bronx_tree_planet_ndvi_time_series.jpg", width = 18, height = 6, units = "in")
