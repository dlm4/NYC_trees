# appeears_downloader_multiple.R

# David Miller, 2024-06-11
# dlm356@cornell.edu

# Modified from:
# https://forum.earthdata.nasa.gov/viewtopic.php?t=2329
# https://github.com/nasa/AppEEARS-Data-Resources/blob/main/guides/How-to-bulk-download-AppEEARS-outputs.md

library(httr)
library(jsonlite)
library(curl)

library(future)
library(future.apply)

# Earthdata username and password
USERNAME <- "your earthdata username"
PASSWORD <- "your earthdata password"

# empty output path for everything, will make directories within this as we go
output_dir <- "/Volumes/NYC_geo/appeears_download/lst_modis_ecostress_20240611/outputs_test/"
# Note: this script will overwrite any files with the same name in a given directory if they exist

# directory only with text files for desired downloads; these are lists from Appeears when select everything and choose to download list
downloadList_dir <- "/Volumes/NYC_geo/appeears_download/lst_modis_ecostress_20240611/download_list/"

#####

# Define the download function
download_file <- function(d, filesl, dataPath) {
  # Create a handle
  s = new_handle()
  handle_setheaders(s, 'Authorization'=paste("Bearer", fromJSON(token_response)$token))
  
  #download
  curl_download(url = filesl[d], destfile = paste0(dataPath, "/", basename(filesl[d])), handle = s)
  Sys.sleep(1)
}

file_list <- list.files(downloadList_dir)

setwd(downloadList_dir)

for (downloadList in file_list){
  print(paste("Working on:", downloadList, sep = " "))
  
  # get list of files from the download list
  filesl <- scan(downloadList, what='list', sep='\n')
  
  setwd(output_dir)
  
  # make new directory and setup
  new_dir_name <- unlist(strsplit(downloadList, "[.]"))[1]
  dir.create(new_dir_name)
  setwd(new_dir_name)
  
  # setup output path
  dataPath <- getwd()
  
  # Create a token by calling AppEEARS API login service. Update the “USERNAME” and “PASSWORD” with yours below
  secret <- base64_enc(paste(USERNAME, PASSWORD, sep = ":"))
  response <- POST("https://appeears.earthdatacloud.nasa.gov/api/login",
                   add_headers("Authorization" = paste("Basic", gsub("\n", "", secret)),
                               "Content-Type" = "application/x-www-form-urlencoded;charset=UTF-8"),
                   body = "grant_type=client_credentials")
  token_response <- prettify(toJSON(content(response), auto_unbox = TRUE))
  
  # Set up the parallel plan
  plan(multisession)  # Use multisession plan for parallel execution
  
  # Execute the download in parallel
  future_lapply(1:length(filesl), function(d) download_file(d, filesl, dataPath))
  
  setwd(downloadList_dir)
}
