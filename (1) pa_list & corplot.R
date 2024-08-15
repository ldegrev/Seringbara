###############################################################################
################################### loading ###################################
###############################################################################

# Set directory
setwd("/Users/lucasdegreve/")

# Load necessary libraries
library(readr)
library(sf)
library(raster)
library(Hmisc)
library(corrplot)

# Load the presence points
nesting <- read_delim("Data/presence/nesting.csv", delim = ",")
feeding_travelling <- read_delim("Data/presence/feeding-travelling.csv", delim = ",")
hunting <- read_delim("Data/presence/hunting.csv", delim = ",")

# Load the tracks and the reserve boundary
tracks <- read_sf("Data/tracks/tracks.shp")
tracks_raster <- raster("Data/tracks_raster.tif")
Nimba_WHS_buffered <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS_buffered.shp")
MCP_buffered <- read_sf("Data/MCP/MCP_buffered.shp")

# Load the distance variables
distance_camp <- raster("GEE/distance/distance_camp.tif")
distance_river <- raster("GEE/distance/distance_river.tif")
distance_road <- raster("GEE/distance/distance_road.tif")
distance_village <- raster("GEE/distance/distance_village.tif")
distance_variables <- c("distance_camp", "distance_river", "distance_road", "distance_village")
for (variable in distance_variables) {
  cropped_raster <- crop(get(variable), extent(MCP_buffered))
  assign(variable, cropped_raster)
}

# Load the environmental variables
raster_directory <- "GEE/Scales"
scales <- c("30m", "100k", "250k", "500k", "1000k", "1500k", "2000k")
variables <- c("aspect", "brightness", "cti", "elev", "evi", "greenness", "hli", "roughness", "slope", "tpi", "treecover", "wetness")
for (scale in scales) {
  for (variable in variables) {
    raster_name <- paste(variable, scale, sep = "_")
    file_path <- file.path(raster_directory, scale, paste(raster_name, "tif", sep = "."))
    raster_layer <- raster(file_path)
    cropped_raster <- crop(raster_layer, extent(MCP_buffered))
    assign(raster_name, cropped_raster)
  }
}

###############################################################################
################################ pseudo-absence ###############################
###############################################################################

################################### feeding ###################################

# Convert tracks_raster to absence points
absence <- rasterToPoints(tracks_raster, spatial = TRUE)
absence_coords <- data.frame(absence@coords)
# Buffer around the presence points
feed_coords = st_as_sf(feeding_travelling, coords=c("X","Y"))
feed_buffer = st_buffer(feed_coords, 30)
absence_sf = st_as_sf(absence_coords, coords=c("x","y"))
b = sapply(st_intersects(absence_sf, feed_buffer), function(x){length(x)==0})
indices0 = which(b, arr.ind = FALSE, useNames = TRUE)
background_pt = absence_coords[c(indices0),]
background_pt_sf = st_as_sf(background_pt,coords=c("x","y"))
# Plotting the buffers
# par(mfrow=c(1,1), mar=c(1,1,1,1), lwd=1)
# plot(distance_road)
# plot(feed_coords, pch = 20, col = "red", add = TRUE)
# plot(feed_buffer, pch = 20, col = rgb(red = 1, green = 0, blue = 0, alpha = 0), add = TRUE)

# Count presence points per raster pixel
coordinates_feed <- data.frame(x = feeding_travelling$X, y = feeding_travelling$Y)
coordinatesSP_feed <- SpatialPoints(coords = coordinates_feed, proj4string = crs(elev_30m))
pointsRaster <- raster(elev_30m) # Initialize a template raster
values(pointsRaster) <- 0 # Set all initial values to 0
pointsCountRaster <- rasterize(coordinatesSP_feed, pointsRaster, fun = 'count')
counts_feed <- as.data.frame(pointsCountRaster, xy = TRUE)
counts_feed <- na.omit(counts_feed)
names(counts_feed) <- c("x", "y", "count")
counts_feed$presence = 1
presence_feed = subset(counts_feed, select = -c(count) ) #Remove column "count"

# Pseudo-absence and weight
weight <- as.data.frame(tracks_raster, xy = TRUE, na.rm = TRUE)
weight$probability <- weight$tracks_raster / sum(weight$tracks_raster) # Calculate probability for each pixel

IDpseudo_df <- data.frame(matrix(ncol = 100, nrow = 479))
pa_feed_list <- list()
# Loop over iterations
for (i in 1:100) {
  IDpseudo <- sample(1:nrow(weight), 479, replace = FALSE, prob = weight$probability)
  IDpseudo_df[, i] <- IDpseudo
  absence_feed <- weight[IDpseudo, 1:2]
  absence_feed$presence <- 0
  pa_feed <- rbind(absence_feed, presence_feed)
  coord_feed <- subset(pa_feed, select = -c(presence))
  pa_feed_list[[i]] <- list(pa_feed = pa_feed, coord_feed = coord_feed)
}
# Save the results list
save(pa_feed_list, file = "pa_feed_list.RData")

################################### hunting ###################################

# Buffer around the presence points
hunt_coords = st_as_sf(hunting, coords=c("X","Y"))
hunt_buffer = st_buffer(hunt_coords, 30)
b = sapply(st_intersects(absence_sf, hunt_buffer), function(x){length(x)==0})
indices0 = which(b, arr.ind = FALSE, useNames = TRUE)
background_pt = absence_coords[c(indices0),]
background_pt_sf = st_as_sf(background_pt,coords=c("x","y"))

# Count presence points per raster pixel
coordinates_hunt <- data.frame(x = hunting$X, y = hunting$Y)
coordinatesSP_hunt <- SpatialPoints(coords = coordinates_hunt, proj4string = crs(elev_30m))
pointsRaster <- raster(elev_30m) # Initialize a template raster
values(pointsRaster) <- 0 # Set all initial values to 0
pointsCountRaster <- rasterize(coordinatesSP_hunt, pointsRaster, fun = 'count')
counts_hunt <- as.data.frame(pointsCountRaster, xy = TRUE)
counts_hunt <- na.omit(counts_hunt)
names(counts_hunt) <- c("x", "y", "count")
counts_hunt$presence = 1
presence_hunt = subset(counts_hunt, select = -c(count) ) #Remove column "count"

IDpseudo_df <- data.frame(matrix(ncol = 100, nrow = 80))
pa_hunt_list <- list()
# Loop over iterations
for (i in 1:100) {
  IDpseudo <- sample(1:nrow(weight), 80, replace = FALSE, prob = weight$probability)
  IDpseudo_df[, i] <- IDpseudo
  absence_hunt <- weight[IDpseudo, 1:2]
  absence_hunt$presence <- 0
  pa_hunt <- rbind(absence_hunt, presence_hunt)
  coord_hunt <- subset(pa_hunt, select = -c(presence))
  pa_hunt_list[[i]] <- list(pa_hunt = pa_hunt, coord_hunt = coord_hunt)
}
# Save the results list
save(pa_hunt_list, file = "pa_hunt_list.RData")

################################### nesting ###################################

# Buffer around the presence points
nest_coords = st_as_sf(nesting, coords=c("X","Y"))
nest_buffer = st_buffer(nest_coords, 30)
b = sapply(st_intersects(absence_sf, nest_buffer), function(x){length(x)==0})
indices0 = which(b, arr.ind = FALSE, useNames = TRUE)
background_pt = absence_coords[c(indices0),]
background_pt_sf = st_as_sf(background_pt,coords=c("x","y"))

# Count presence points per raster pixel
coordinates_nest <- data.frame(x = nesting$X, y = nesting$Y)
coordinatesSP_nest <- SpatialPoints(coords = coordinates_nest, proj4string = crs(elev_30m))
pointsRaster <- raster(elev_30m) # Initialize a template raster
values(pointsRaster) <- 0 # Set all initial values to 0
pointsCountRaster <- rasterize(coordinatesSP_nest, pointsRaster, fun = 'count')
counts_nest <- as.data.frame(pointsCountRaster, xy = TRUE)
counts_nest <- na.omit(counts_nest)
names(counts_nest) <- c("x", "y", "count")
counts_nest$presence = 1
presence_nest = subset(counts_nest, select = -c(count) ) #Remove column "count"

IDpseudo_df <- data.frame(matrix(ncol = 100, nrow = 174))
pa_nest_list <- list()
# Loop over iterations
for (i in 1:100) {
  IDpseudo <- sample(1:nrow(weight), 174, replace = FALSE, prob = weight$probability)
  IDpseudo_df[, i] <- IDpseudo
  absence_nest <- weight[IDpseudo, 1:2]
  absence_nest$presence <- 0
  pa_nest <- rbind(absence_nest, presence_nest)
  coord_nest <- subset(pa_nest, select = -c(presence))
  pa_nest_list[[i]] <- list(pa_nest = pa_nest, coord_nest = coord_nest)
}
# Save the results list
save(pa_nest_list, file = "pa_nest_list.RData")


###############################################################################
################################# correlation #################################
###############################################################################

################################### feeding ###################################

library(Hmisc)
library(corrplot)

load("Data/pa_list/pa_feed_list.RData")

cor_list <- list()
for (i in 1:100) { # Loop through the 100th coord_feed elements
  extracted_data_feed <- list()
  for (variable in variables) {   # Extract data for variables
    raster_layer <- get(paste0(variable, "_30m"))
    extracted_data_feed[[variable]] <- extract(raster_layer, pa_feed_list[[i]][["coord_feed"]])
  }
  for (variable in distance_variables) {   # Extract data for distance variables
    raster_layer <- get(variable)
    extracted_data_feed[[variable]] <- extract(raster_layer, pa_feed_list[[i]][["coord_feed"]])
  }
  data_feed <- cbind(do.call(cbind, extracted_data_feed))   # Combine all extracted data into a dataframe
  env.cor <- cor(data_feed, method = "spearman")   # Compute the correlation matrix using Spearman method
  cor_list[[i]] <- env.cor   # Store the correlation matrix in the list
}
avg_cor_matrix_feed <- Reduce("+", cor_list) / length(cor_list) # Compute the average correlation matrix
save(avg_cor_matrix_feed, file = "avg_cor_matrix_feed.RData")
corrplot(avg_cor_matrix_feed, method = "color", type = "upper", addCoef.col = 1, number.cex = 0.65, pch.cex = 1, diag = FALSE) # Plot the average correlation matrix

# load("/Users/lucasdegreve/Data/RData/cor_matrix/avg_cor_matrix_feed.RData")
# cor.score <- rowSums(abs(avg_cor_matrix_feed))
# cor.score.ranked <- sort(cor.score, decreasing = TRUE)
# print(cor.score.ranked)

################################### hunting ###################################

load("Data/pa_list/pa_hunt_list.RData")

cor_list <- list()
for (i in 1:100) { # Loop through the 100th coord_hunt elements
  extracted_data_hunt <- list()
  for (variable in variables) {   # Extract data for variables
    raster_layer <- get(paste0(variable, "_30m"))
    extracted_data_hunt[[variable]] <- extract(raster_layer, pa_hunt_list[[i]][["coord_hunt"]])
  }
  for (variable in distance_variables) {   # Extract data for distance variables
    raster_layer <- get(variable)
    extracted_data_hunt[[variable]] <- extract(raster_layer, pa_hunt_list[[i]][["coord_hunt"]])
  }
  data_hunt <- cbind(do.call(cbind, extracted_data_hunt))   # Combine all extracted data into a dataframe
  env.cor <- cor(data_hunt, method = "spearman")   # Compute the correlation matrix using Spearman method
  cor_list[[i]] <- env.cor   # Store the correlation matrix in the list
}
avg_cor_matrix_hunt <- Reduce("+", cor_list) / length(cor_list) # Compute the average correlation matrix
save(avg_cor_matrix_hunt, file = "avg_cor_matrix_hunt.RData")
corrplot(avg_cor_matrix_hunt, method = "color", type = "upper", addCoef.col = 1, number.cex = 0.65, pch.cex = 1, diag = FALSE) # Plot the average correlation matrix

################################### nesting ###################################

load("Data/pa_list/pa_nest_list.RData")

cor_list <- list()
for (i in 1:100) { # Loop through the 100th coord_nest elements
  extracted_data_nest <- list()
  for (variable in variables) {   # Extract data for variables
    raster_layer <- get(paste0(variable, "_30m"))
    extracted_data_nest[[variable]] <- extract(raster_layer, pa_nest_list[[i]][["coord_nest"]])
  }
  for (variable in distance_variables) {   # Extract data for distance variables
    raster_layer <- get(variable)
    extracted_data_nest[[variable]] <- extract(raster_layer, pa_nest_list[[i]][["coord_nest"]])
  }
  data_nest <- cbind(do.call(cbind, extracted_data_nest))   # Combine all extracted data into a dataframe
  env.cor <- cor(data_nest, method = "spearman")   # Compute the correlation matrix using Spearman method
  cor_list[[i]] <- env.cor   # Store the correlation matrix in the list
}
avg_cor_matrix_nest <- Reduce("+", cor_list) / length(cor_list) # Compute the average correlation matrix
save(avg_cor_matrix_nest, file = "avg_cor_matrix_nest.RData")
corrplot(avg_cor_matrix_nest, method = "color", type = "upper", addCoef.col = 1, number.cex = 0.65, pch.cex = 1, diag = FALSE) # Plot the average correlation matrix

################################## centroids ##################################

# Create centroids
centroids_MCP <- rasterToPoints(elev_30m, spatial = TRUE)
centroids_MCP_df <- data.frame(centroids_MCP@coords)

extracted_data_MCP <- list()
for (variable in variables) { # Extract data for variables
  raster_layer <- get(paste0(variable, "_30m"))
  extracted_data_MCP[[variable]] <- extract(raster_layer, centroids_MCP_df)
}
for (variable in distance_variables) { # Extract data for distance variables
  raster_layer <- get(variable)
  extracted_data_MCP[[variable]] <- extract(raster_layer, centroids_MCP_df)
}
# Combine all extracted data into a dataframe
data_raster_MCP <- cbind(do.call(cbind, extracted_data_MCP))

env.cor = cor(data_raster_MCP, method = c("spearman"))
env.rcorr = rcorr(as.matrix(data_raster_MCP))
env.coeff = env.rcorr$r #r : the correlation matrix
env.p = env.rcorr$P #P : the p-values corresponding to the significance levels of correlations
corrplot(env.cor, method = "color", type = "upper", addCoef.col = 1, number.cex = 0.65, pch.cex = 1, diag = FALSE) # Plot the average correlation matrix
save(env.cor, file = "cor_matrix_MCP.RData")

###############################################################################
################################## 85% lists ##################################
###############################################################################

library(dplyr)
# Function to split data into 85% and 15% sets while maintaining class balance
split_data <- function(pa_df, coord_df) {
  # Split presence and absence separately
  pres <- pa_df %>% filter(presence == 1)
  abs <- pa_df %>% filter(presence == 0)
  # Determine the split indices
  pres_split <- sample(nrow(pres), size = 0.85 * nrow(pres))
  abs_split <- sample(nrow(abs), size = 0.85 * nrow(abs))
  # Create 85% and 15% splits
  pa_85 <- bind_rows(pres[pres_split, ], abs[abs_split, ])
  pa_15 <- bind_rows(pres[-pres_split, ], abs[-abs_split, ])
  # Corresponding coord splits
  coord_85 <- coord_df[row.names(pa_85), ]
  coord_15 <- coord_df[row.names(pa_15), ]
  return(list(pa_85 = pa_85, coord_85 = coord_85, pa_15 = pa_15, coord_15 = coord_15))
}

# Function to process a list
process_list <- function(pa_list, type) {
  new_list <- list()
  for (i in 1:length(pa_list)) {
    pa_data <- pa_list[[i]][[paste0("pa_", type)]]
    coord_data <- pa_list[[i]][[paste0("coord_", type)]]
    split_result <- split_data(pa_data, coord_data)
    new_list[[i]] <- list()
    new_list[[i]][[paste0("pa_", type, "85")]] <- split_result$pa_85
    new_list[[i]][[paste0("coord_", type, "85")]] <- split_result$coord_85
    new_list[[i]][[paste0("pa_", type, "15")]] <- split_result$pa_15
    new_list[[i]][[paste0("coord_", type, "15")]] <- split_result$coord_15
  }
  return(new_list)
}

# Load the data
load("/Users/lucasdegreve/Data/RData/pa_list/pa_feed_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_hunt_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_nest_list.RData")

# Process each list
pa_feed_list85 <- process_list(pa_feed_list, "feed")
pa_hunt_list85 <- process_list(pa_hunt_list, "hunt")
pa_nest_list85 <- process_list(pa_nest_list, "nest")

# Save the new lists if needed
save(pa_feed_list85, file = "pa_feed_list85.RData")
save(pa_hunt_list85, file = "pa_hunt_list85.RData")
save(pa_nest_list85, file = "pa_nest_list85.RData")