
###############################################################################
################################### loading ###################################
###############################################################################

# Set directory
setwd("/Users/lucasdegreve/")

# Load necessary libraries
library(raster)
library(sf)
library(MuMIn)
library(fields)
library(RColorBrewer)
library(dplyr)
library(openxlsx)

# Load Nimba boudary and centroids_df
Nimba_WHS <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS.shp")
Nimba_WHS_buffered <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS_buffered.shp")
MCP_buffered <- read_sf("Data/MCP/MCP_buffered.shp")
MCP <- read_sf("Data/MCP/MCP.shp")
elev_30m <- raster("GEE/30m/elev_30m.tif")
elev_30m <- crop(elev_30m, extent(Nimba_WHS_buffered))
centroids <- rasterToPoints(elev_30m, spatial = TRUE)
centroids_df <- data.frame(centroids@coords)

# Load the distance variables
distance_camp <- raster("GEE/distance/distance_camp.tif")
distance_river <- raster("GEE/distance/distance_river.tif")
distance_road <- raster("GEE/distance/distance_road.tif")
distance_village <- raster("GEE/distance/distance_village.tif")
distance_variables <- c("distance_camp", "distance_river", "distance_road", "distance_village")
for (variable in distance_variables) {
  cropped_raster <- crop(get(variable), extent(Nimba_WHS_buffered))
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
    cropped_raster <- crop(raster_layer, extent(Nimba_WHS_buffered))
    assign(raster_name, cropped_raster)
  }
}

###############################################################################
############################### GLM data_raster ############################### 
###############################################################################

# Define the function to create data frames for each GLM model
create_data_raster <- function(variable_names, rasters, centroids_df) {
  extracted_data <- lapply(rasters, extract, y = centroids_df)
  data_raster <- cbind(centroids_df, do.call(cbind, extracted_data))
  colnames(data_raster) <- c("x", "y", variable_names)
  return(data_raster)
}

# Define the variable names for each GLM model type
glm_variable_names <- list(
  glm_feed_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  glm_hunt_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  glm_nest_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  glm_feed_scale = c("aspect_2000k", "cti_1000k", "elev_250k", "evi_1500k", "hli_1000k", "slope_250k", "tpi_2000k", "treecover_1500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  glm_hunt_scale = c("aspect_1000k", "cti_500k", "elev_30m", "evi_1500k", "slope_1000k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  glm_nest_scale = c("aspect_1500k", "cti_1000k", "elev_250k", "evi_1500k", "hli_100k", "slope_250k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
)

# Define the lists of rasters for each GLM model type
glm_rasters <- list(
  glm_feed_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  glm_hunt_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  glm_nest_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  glm_feed_scale = list(aspect_2000k, cti_1000k, elev_250k, evi_1500k, hli_1000k, slope_250k, tpi_2000k, treecover_1500k, wetness_1500k, distance_camp, distance_river, distance_road, distance_village),
  glm_hunt_scale = list(aspect_1000k, cti_500k, elev_30m, evi_1500k, slope_1000k, tpi_2000k, treecover_1000k, wetness_1500k, distance_camp, distance_river, distance_road, distance_village),
  glm_nest_scale = list(aspect_1500k, cti_1000k, elev_250k, evi_1500k, hli_100k, slope_250k, tpi_2000k, treecover_1000k, wetness_1500k, distance_camp, distance_river, distance_road, distance_village)
)

# Create the data frames for each GLM model
data_rasters <- lapply(names(glm_variable_names), function(name) {
  create_data_raster(glm_variable_names[[name]], glm_rasters[[name]], centroids_df)
})

# Assign data frames to variables
data_raster_glm_feed_30m <- data_rasters[[1]]
data_raster_glm_hunt_30m <- data_rasters[[2]]
data_raster_glm_nest_30m <- data_rasters[[3]]
data_raster_glm_feed_scale <- data_rasters[[4]]
data_raster_glm_hunt_scale <- data_rasters[[5]]
data_raster_glm_nest_scale <- data_rasters[[6]]

save(data_raster_glm_feed_30m, file = "data_raster_glm_feed_30m.RData")
save(data_raster_glm_hunt_30m, file = "data_raster_glm_hunt_30m.RData")
save(data_raster_glm_nest_30m, file = "data_raster_glm_nest_30m.RData")
save(data_raster_glm_feed_scale, file = "data_raster_glm_feed_scale.RData")
save(data_raster_glm_hunt_scale, file = "data_raster_glm_hunt_scale.RData")
save(data_raster_glm_nest_scale, file = "data_raster_glm_nest_scale.RData")

###############################################################################
############################### GLM prediction ################################
###############################################################################

load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_nest_scale.RData")

load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_feed_30m_list85_dredge_results.RData")
load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_feed_scale_list85_dredge_results.RData")
load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_hunt_30m_list85_dredge_results.RData")
load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_hunt_scale_list85_dredge_results.RData")
load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_nest_30m_list85_dredge_results.RData")
load("/Users/lucasdegreve/Data/RData/glm_results/df_glm_nest_scale_list85_dredge_results.RData")

load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_feed_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_feed_scale_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_hunt_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_hunt_scale_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_nest_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_nest_scale_list85.RData")

avg_list_feed_30m <- list()
predictions_df_glm_feed_30m <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_feed_30m_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_feed_30m_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = FALSE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_feed_30m)
  avg_list_feed_30m[[i]] <- avg
  predictions_df_glm_feed_30m[, i] <- predictions
}
save(avg_list_feed_30m, file = "avg_list_feed_30m.RData")
save(predictions_df_glm_feed_30m, file = "predictions_df_glm_feed_30m.RData")

avg_list_feed_scale <- list()
predictions_df_glm_feed_scale <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_feed_scale_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_feed_scale_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = FALSE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_feed_scale)
  avg_list_feed_scale[[i]] <- avg
  predictions_df_glm_feed_scale[, i] <- predictions
}
save(avg_list_feed_scale, file = "avg_list_feed_scale.RData")
save(predictions_df_glm_feed_scale, file = "predictions_df_glm_feed_scale.RData")

avg_list_hunt_30m <- list()
predictions_df_glm_hunt_30m <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_hunt_30m_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_hunt_30m_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = TRUE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_hunt_30m)
  avg_list_hunt_30m[[i]] <- avg
  predictions_df_glm_hunt_30m[, i] <- predictions
}
save(avg_list_hunt_30m, file = "avg_list_hunt_30m.RData")
save(predictions_df_glm_hunt_30m, file = "predictions_df_glm_hunt_30m.RData")

avg_list_hunt_scale <- list()
predictions_df_glm_hunt_scale <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_hunt_scale_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_hunt_scale_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = FALSE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_hunt_scale)
  avg_list_hunt_scale[[i]] <- avg
  predictions_df_glm_hunt_scale[, i] <- predictions
}
save(avg_list_hunt_scale, file = "avg_list_hunt_scale.RData")
save(predictions_df_glm_hunt_scale, file = "predictions_df_glm_hunt_scale.RData")

avg_list_nest_30m <- list()
predictions_df_glm_nest_30m <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_nest_30m_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_nest_30m_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = FALSE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_nest_30m)
  avg_list_nest_30m[[i]] <- avg
  predictions_df_glm_nest_30m[, i] <- predictions
}
save(avg_list_nest_30m, file = "avg_list_nest_30m.RData")
save(predictions_df_glm_nest_30m, file = "predictions_df_glm_nest_30m.RData")

avg_list_nest_scale <- list()
predictions_df_glm_nest_scale <- data.frame(matrix(ncol = 100, nrow = 666490))
for (i in 1:100) {
  df <- df_glm_nest_scale_list85[[i]]
  formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
  global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
  # Get the dredge results for the current model
  dredge_results_dup <- df_glm_nest_scale_list85_dredge_results[[i]]
  # Filter models with delta < 2
  models_under_delta2 <- dredge_results_dup[dredge_results_dup$delta < 2, ]
  models_under_delta2_names <- rownames(models_under_delta2)
  # Create the dredge results list
  dredge_results_list <- lapply(dredge(global_model, trace = FALSE, evaluate = FALSE), eval)
  dredge_results_list_filtered <- dredge_results_list[names(dredge_results_list) %in% models_under_delta2_names]
  avg <- model.avg(dredge_results_list_filtered)
  predictions <- predict(avg, newdata = data_raster_glm_nest_scale)
  avg_list_nest_scale[[i]] <- avg
  predictions_df_glm_nest_scale[, i] <- predictions
}
save(avg_list_nest_scale, file = "avg_list_nest_scale.RData")
save(predictions_df_glm_nest_scale, file = "predictions_df_glm_nest_scale.RData")

###############################################################################
################################### GLM maps ################################## 
###############################################################################

# Load necessary libraries
library(RColorBrewer)
library(fields)
library(sf)
library(dplyr)

# Load the MCP boundary
elev_30m_MCP <- crop(elev_30m, extent(MCP_buffered))
load("/Users/lucasdegreve/Data/RData/centroids_plot.RData")
centroids_sf <- st_as_sf(centroids_df, coords = c("x", "y"), crs = st_crs(Nimba_WHS_buffered))
MCP_buffered <- st_transform(MCP_buffered, st_crs(centroids_sf))
inside_MCP <- st_within(centroids_sf, MCP_buffered, sparse = FALSE)
centroids_MCP <- centroids_plot
centroids_MCP[!inside_MCP, c("x", "y")] <- NA
Nimba_WHS_within_MCP <- st_intersection(Nimba_WHS, MCP_buffered)
valid_indices <- !is.na(centroids_MCP$x) & !is.na(centroids_MCP$y)
filtered_centroids_MCP <- centroids_MCP[valid_indices, ]

# Define plot dimensions
plot_extent <- extent(elev_30m_MCP)
plot_width <- plot_extent@xmax - plot_extent@xmin
plot_height <- plot_extent@ymax - plot_extent@ymin

# Define the plotting function
plot_predictions <- function(predictions_df, centroids_MCP, valid_indices, model_name, plot_width, plot_height, Nimba_WHS_within_MCP, MCP_buffered) {
  # Define output file names
  mean_plot_file <- paste0("plot_", model_name, "_mean.png")
  sd_plot_file <- paste0("plot_", model_name, "_sd.png")
  # Plot mean predictions
  png(mean_plot_file, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  prediction_means_prob <- 1 / (1 + exp(-predictions_df))
  prediction_means <- rowMeans(prediction_means_prob)
  filtered_prediction_means <- prediction_means[valid_indices]
  min_means <- min(filtered_prediction_means)
  max_means <- max(filtered_prediction_means)
  colourScale = colorRampPalette(brewer.pal(11, "YlGn"))(121)[11:121]
  cols <- colourScale[(((prediction_means - min_means) / (max_means - min_means)) * 100) + 1]
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", main = paste("Prediction -", model_name), frame.plot = FALSE)
  image.plot(z = filtered_prediction_means, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
  
  # Plot standard deviation
  png(sd_plot_file, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  prediction_sd <- apply(predictions_df, 1, sd)
  prediction_sd_log <- 1 / (1 + exp(-prediction_sd))
  prediction_sd1 <- 1 / (1 + exp(-predictions_df))
  prediction_sd1_log <- apply(prediction_sd1, 1, sd)
  filtered_prediction_sd <- prediction_sd_log[valid_indices]
  min_sd <- min(filtered_prediction_sd)
  max_sd <- max(filtered_prediction_sd)
  colourScale <- colorRampPalette(brewer.pal(11, "OrRd"))(121)[11:121]
  cols <- colourScale[(((prediction_sd_log - min_sd) / (max_sd - min_sd)) * 100) + 1]
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", main = paste("Standard deviation -", model_name), frame.plot = FALSE)
  image.plot(z = prediction_sd1_log, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
}

# List of models and their corresponding names
model_files <- list(
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_feed_30m.RData", name = "feed_30m"),
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_feed_scale.RData", name = "feed_scale"),
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_hunt_30m.RData", name = "hunt_30m"),
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_hunt_scale.RData", name = "hunt_scale"),
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_nest_30m.RData", name = "nest_30m"),
  list(file = "/Users/lucasdegreve/Data/RData/prediction_df/glm/predictions_df_glm_nest_scale.RData", name = "nest_scale")
)

# Loop through the models
for (model_info in model_files) {
  # Load the model data
  load(model_info$file)
  # Extract the predictions dataframe variable
  predictions_df <- get(ls(pattern = paste0("predictions_df_glm_", model_info$name)))
  # Apply the plotting function
  plot_predictions(predictions_df, centroids_MCP, valid_indices, model_info$name, plot_width, plot_height, Nimba_WHS_within_MCP, MCP_buffered)
}

###############################################################################
################################### Response ################################## 
###############################################################################

load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/data_raster/glm/data_raster_glm_nest_scale.RData")

load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_scale.RData")

# Load the MCP boundary
elev_30m_MCP <- crop(elev_30m, extent(MCP_buffered))
load("/Users/lucasdegreve/Data/RData/centroids_plot.RData")
centroids_sf <- st_as_sf(centroids_df, coords = c("x", "y"), crs = st_crs(Nimba_WHS_buffered))
MCP_buffered <- st_transform(MCP_buffered, st_crs(centroids_sf))
inside_MCP <- st_within(centroids_sf, MCP_buffered, sparse = FALSE)
centroids_MCP <- centroids_plot
centroids_MCP[!inside_MCP, c("x", "y")] <- NA
Nimba_WHS_within_MCP <- st_intersection(Nimba_WHS, MCP_buffered)

# Define plot dimensions
plot_extent <- extent(elev_30m_MCP)
plot_width <- plot_extent@xmax - plot_extent@xmin
plot_height <- plot_extent@ymax - plot_extent@ymin

# Define a function to process each model
process_model <- function(model_name) {
  # Load data
  data_raster <- get(paste0("data_raster_glm_", model_name))
  avg_list <- get(paste0("avg_list_", model_name))
  # Initialize predictions dataframe
  predictions_df <- data.frame(matrix(ncol = 100, nrow = nrow(data_raster)))
  # Generate predictions
  for (i in 1:100) {
    avg <- avg_list[[i]]
    predictions <- predict(avg, newdata = data_raster, type = "response")
    predictions_df[, i] <- predictions
  }
  # Save predictions dataframe
  save(predictions_df, file = paste0("predictions_df_glm_", model_name, ".RData"))
  # Calculate means
  prediction_means <- rowMeans(predictions_df)
  min_means <- min(prediction_means)
  max_means <- max(prediction_means)
  colourScale = colorRampPalette(brewer.pal(11, "YlGn"))(121)[11:121]
  cols <- colourScale[(((prediction_means - min_means) / (max_means - min_means)) * 100) + 1]
  # Plot means
  png(paste0(model_name, "_means.png"), width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", frame.plot = FALSE)
  image.plot(z = prediction_means, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()
  # Calculate standard deviation
  prediction_sd <- apply(predictions_df, 1, sd)
  min_sd <- min(prediction_sd)
  max_sd <- max(prediction_sd)
  colourScale <- colorRampPalette(brewer.pal(11, "OrRd"))(121)[11:121]
  cols <- colourScale[(((prediction_sd - min_sd) / (max_sd - min_sd)) * 100) + 1]
  # Plot standard deviation
  png(paste0(model_name, "_sd.png"), width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", frame.plot = FALSE)
  image.plot(z = prediction_sd, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()
}

# Apply the function to each model
model_names <- c("feed_30m", "feed_scale", "hunt_30m", "hunt_scale", "nest_30m", "nest_scale")
lapply(model_names, process_model)

###############################################################################
################################ 15% Dataframes ############################### 
###############################################################################

# Load the MCP boundary
MCP_buffered <- read_sf("Data/MCP/MCP_buffered.shp")
Nimba_WHS <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS.shp")
Nimba_WHS_buffered <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS_buffered.shp")

# Load the distance variables
distance_camp <- raster("GEE/distance/distance_camp.tif")
distance_river <- raster("GEE/distance/distance_river.tif")
distance_road <- raster("GEE/distance/distance_road.tif")
distance_village <- raster("GEE/distance/distance_village.tif")
distance_variables <- c("distance_camp", "distance_river", "distance_road", "distance_village")
for (variable in distance_variables) {
  cropped_raster <- crop(get(variable), extent(Nimba_WHS_buffered))
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
    cropped_raster <- crop(raster_layer, extent(Nimba_WHS_buffered))
    assign(raster_name, cropped_raster)
  }
}

# Function to create data frames for a given list, type, and variable set
create_data_frames_v2 <- function(pa_list, type, variables) {
  df_list <- list()
  for (i in 1:100) {
    pa_data <- pa_list[[i]][[paste0("pa_", type)]]
    coord_data <- pa_list[[i]][[paste0("coord_", type)]]
    scale_list <- list()
    variable_names <- c("x", "y", "presence")
    for (variable in variables) {
      raster_data <- extract(get(variable), coord_data)
      scale_list[[length(scale_list) + 1]] <- raster_data
      variable_names <- c(variable_names, variable)
    }
    df <- cbind(pa_data, do.call(cbind, scale_list))
    colnames(df) <- variable_names
    df_list[[i]] <- df
  }
  return(df_list)
}

load("/Users/lucasdegreve/Data/RData/pa_list/pa_hunt_list85.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_nest_list85.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_feed_list85.RData")

# correlation filter 0.75 / 0.8
glm_feed_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_hunt_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_nest_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")

# correlation filter 0.77 / 0.8 (-hli_1500k for hunt model)
glm_feed_scale <- c("aspect_2000k", "cti_1000k", "elev_250k", "evi_1500k", "hli_1000k", "slope_250k", "tpi_2000k", "treecover_1500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_hunt_scale <- c("aspect_1000k", "cti_500k", "elev_30m", "evi_1500k", "slope_1000k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_nest_scale <- c("aspect_1500k", "cti_1000k", "elev_250k", "evi_1500k", "hli_100k", "slope_250k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")

df_glm_feed_30m_list15 <- create_data_frames_v2 (pa_feed_list85, "feed15", glm_feed_30m)
df_glm_hunt_30m_list15 <- create_data_frames_v2 (pa_hunt_list85, "hunt15", glm_hunt_30m)
df_glm_nest_30m_list15 <- create_data_frames_v2 (pa_nest_list85, "nest15", glm_nest_30m)

df_glm_feed_scale_list15 <- create_data_frames_v2 (pa_feed_list85, "feed15", glm_feed_scale)
df_glm_hunt_scale_list15 <- create_data_frames_v2 (pa_hunt_list85, "hunt15", glm_hunt_scale)
df_glm_nest_scale_list15 <- create_data_frames_v2 (pa_nest_list85, "nest15", glm_nest_scale)

df_15_lists <- list(
  df_glm_feed_30m_list15 = df_glm_feed_30m_list15,
  df_glm_hunt_30m_list15 = df_glm_hunt_30m_list15,
  df_glm_nest_30m_list15 = df_glm_nest_30m_list15,
  df_glm_feed_scale_list15 = df_glm_feed_scale_list15,
  df_glm_hunt_scale_list15 = df_glm_hunt_scale_list15,
  df_glm_nest_scale_list15 = df_glm_nest_scale_list15
)

# Loop over the list and save each item
for (name in names(df_15_lists)) {
  save(list = name, file = paste0(name, ".RData"))
}

###############################################################################
########################### GLM Results & Validation ########################## 
###############################################################################

library(pROC)
library(openxlsx)

# Load the datasets
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_feed_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_feed_scale_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_hunt_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_hunt_scale_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_nest_30m_list85.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/85/df_glm_nest_scale_list85.RData")

load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_feed_30m_list15.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_feed_scale_list15.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_hunt_30m_list15.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_hunt_scale_list15.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_nest_30m_list15.RData")
load("/Users/lucasdegreve/Data/RData/df_glm/15/df_glm_nest_scale_list15.RData")

load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_scale.RData")

summary(avg_list_hunt_scale[[1]])

# Define the models
models <- list(
  feed_30m = list(avg_list = avg_list_feed_30m, df_list85 = df_glm_feed_30m_list85, df_list15 = df_glm_feed_30m_list15),
  feed_scale = list(avg_list = avg_list_feed_scale, df_list85 = df_glm_feed_scale_list85, df_list15 = df_glm_feed_scale_list15),
  hunt_30m = list(avg_list = avg_list_hunt_30m, df_list85 = df_glm_hunt_30m_list85, df_list15 = df_glm_hunt_30m_list15),
  hunt_scale = list(avg_list = avg_list_hunt_scale, df_list85 = df_glm_hunt_scale_list85, df_list15 = df_glm_hunt_scale_list15),
  nest_30m = list(avg_list = avg_list_nest_30m, df_list85 = df_glm_nest_30m_list85, df_list15 = df_glm_nest_30m_list15),
  nest_scale = list(avg_list = avg_list_nest_scale, df_list85 = df_glm_nest_scale_list85, df_list15 = df_glm_nest_scale_list15)
)

# Process each model
process_model <- function(model_name, model_data) {
  aic_values <- numeric(100)
  auc_values_train <- numeric(100)
  auc_values_test <- numeric(100)
  r_squared_train <- numeric(100)
  r_squared_test <- numeric(100)
  results_list <- list()
  
  for (i in 1:100) {
    avg <- model_data$avg_list[[i]]
    temp_results_df <- summary(avg)$coefmat.subset
    results_list[[i]] <- temp_results_df
    
    aic_value <- mean(summary(avg)[["msTable"]][["AICc"]])
    aic_values[i] <- aic_value
    
    predictions_train <- predict(avg, newdata = model_data$df_list85[[i]], type = "response")
    auc_value_train <- roc(model_data$df_list85[[i]]$presence, predictions_train)$auc
    auc_values_train[i] <- auc_value_train
    
    observed_train <- model_data$df_list85[[i]]$presence
    residuals_train <- observed_train - predictions_train
    ss_res_train <- sum(residuals_train^2)
    ss_tot_train <- sum((observed_train - mean(observed_train))^2)
    r_squared_train[i] <- 1 - (ss_res_train / ss_tot_train)
    
    predictions_test <- predict(avg, newdata = model_data$df_list15[[i]], type = "response")
    auc_value_test <- roc(model_data$df_list15[[i]]$presence, predictions_test)$auc
    auc_values_test[i] <- auc_value_test
    
    observed_test <- model_data$df_list15[[i]]$presence
    residuals_test <- observed_test - predictions_test
    ss_res_test <- sum(residuals_test^2)
    ss_tot_test <- sum((observed_test - mean(observed_test))^2)
    r_squared_test[i] <- 1 - (ss_res_test / ss_tot_test)
  }
  
  all_vars_stats <- list()
  for (i in 1:length(results_list)) {
    temp_results_df <- results_list[[i]]
    colnames(temp_results_df) <- c("Estimate", "Std. Error", "Adjusted SE", "z value", "Pr(>|z|)")
    for (var in rownames(temp_results_df)) {
      if (!(var %in% names(all_vars_stats))) {
        all_vars_stats[[var]] <- data.frame(Estimate = numeric(),
                                            `Std. Error` = numeric(),
                                            `Adjusted SE` = numeric(),
                                            `z value` = numeric(),
                                            `Pr(>|z|)` = numeric())
      }
      all_vars_stats[[var]] <- rbind(all_vars_stats[[var]], temp_results_df[var, ])
      colnames(all_vars_stats[[var]]) <- c("Estimate", "Std. Error", "Adjusted SE", "z value", "Pr(>|z|)")
    }
  }
  
  calc_means <- function(df) {
    data.frame(Estimate = mean(df$Estimate, na.rm = TRUE),
               `Std. Error` = mean(df$`Std. Error`, na.rm = TRUE),
               `Adjusted SE` = mean(df$`Adjusted SE`, na.rm = TRUE),
               `z value` = mean(df$`z value`, na.rm = TRUE),
               `Pr(>|z|)` = mean(df$`Pr(>|z|)`, na.rm = TRUE))
  }
  mean_stats_list <- lapply(all_vars_stats, calc_means)
  results_df <- do.call(rbind, mean_stats_list)
  results_df <- data.frame(Variable = rownames(results_df), results_df, row.names = NULL)
  
  write.xlsx(results_df, file = paste0("results_df_", model_name, ".xlsx"), sheetName = "Results GLM", colNames = TRUE)
  
  summary_stats <- data.frame(
    Metric = c("AIC", "AUC Train", "R2 Train", "AUC Test", "R2 Test"),
    Min = c(min(aic_values), min(auc_values_train), min(r_squared_train), min(auc_values_test), min(r_squared_test)),
    Max = c(max(aic_values), max(auc_values_train), max(r_squared_train), max(auc_values_test), max(r_squared_test)),
    Mean = c(mean(aic_values), mean(auc_values_train), mean(r_squared_train), mean(auc_values_test), mean(r_squared_test)),
    SD = c(sd(aic_values), sd(auc_values_train), sd(r_squared_train), sd(auc_values_test), sd(r_squared_test))
  )
  
  write.xlsx(summary_stats, file = paste0("summary_stats_", model_name, ".xlsx"), sheetName = "Summary Stats", colNames = TRUE)
}

# Loop through all models
for (model_name in names(models)) {
  process_model(model_name, models[[model_name]])
}

###############################################################################
################################ BRT R squared ################################
###############################################################################

load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_scale_list_brt_results.RData")

load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_feed_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_feed_scale_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_hunt_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_hunt_scale_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_nest_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_nest_scale_list.RData")

library(gbm)
library(openxlsx)
# Function to calculate R-squared for a list of BRT models and corresponding data frames
calculate_r_squared <- function(brt_models, data_list) {
  r_squared_values <- numeric(length(brt_models))
  for (i in 1:length(brt_models)) {
    brt_model <- brt_models[[i]]
    prediction <- predict.gbm(brt_model, data_list[[i]], brt_model$gbm.call$best.trees, type = "response", single.tree = FALSE)
    observed <- data_list[[i]]$presence
    mean_observed <- mean(observed)
    sst <- sum((observed - mean_observed)^2)
    ssr <- sum((observed - prediction)^2)
    r_squared <- 1 - (ssr / sst)
    r_squared_values[i] <- r_squared
  }
  return(r_squared_values)
}
# Load BRT models and data lists
brt_results_list <- list(
  feed_30m = df_brt_feed_30m_list_brt_results,
  feed_scale = df_brt_feed_scale_list_brt_results,
  hunt_30m = df_brt_hunt_30m_list_brt_results,
  hunt_scale = df_brt_hunt_scale_list_brt_results,
  nest_30m = df_brt_nest_30m_list_brt_results,
  nest_scale = df_brt_nest_scale_list_brt_results
)
data_lists <- list(
  feed_30m = df_brt_feed_30m_list,
  feed_scale = df_brt_feed_scale_list,
  hunt_30m = df_brt_hunt_30m_list,
  hunt_scale = df_brt_hunt_scale_list,
  nest_30m = df_brt_nest_30m_list,
  nest_scale = df_brt_nest_scale_list
)
# Initialize a list to store summary statistics for each model
summary_stats_list <- list()
# Loop over each model to calculate R-squared and summary statistics
for (model_name in names(brt_results_list)) {
  r_squared_values <- calculate_r_squared(brt_results_list[[model_name]], data_lists[[model_name]])
  summary_stats <- data.frame(
    Metric = "R-squared",
    Min = min(r_squared_values, na.rm = TRUE),
    Max = max(r_squared_values, na.rm = TRUE),
    Mean = mean(r_squared_values, na.rm = TRUE),
    SD = sd(r_squared_values, na.rm = TRUE)
  )
  summary_stats_list[[model_name]] <- summary_stats
}
# Combine all summary statistics into a single data frame
final_summary_stats <- do.call(rbind, summary_stats_list)
rownames(final_summary_stats) <- names(summary_stats_list)
write.xlsx(final_summary_stats, file = "summary_stats_brt_models.xlsx", sheetName = "Summary Stats BRT Models", colNames = TRUE)

###############################################################################
################################ BRT train AUC ################################
###############################################################################

# Load the required datasets
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_scale_list_brt_results.RData")

# List of models
brt_model_lists <- list(
  feed_30m = df_brt_feed_30m_list_brt_results,
  feed_scale = df_brt_feed_scale_list_brt_results,
  hunt_30m = df_brt_hunt_30m_list_brt_results,
  hunt_scale = df_brt_hunt_scale_list_brt_results,
  nest_30m = df_brt_nest_30m_list_brt_results,
  nest_scale = df_brt_nest_scale_list_brt_results
)

# Initialize a vector to store the average AUCs
avg_aucs <- numeric(length(brt_model_lists))
names(avg_aucs) <- names(brt_model_lists)

# Calculate the average AUC for each model list
for (model_name in names(brt_model_lists)) {
  auc_values_train <- numeric(100)
  for (i in 1:100) {
    auc_value_train <- brt_model_lists[[model_name]][[i]][["self.statistics"]][["discrimination"]]
    auc_values_train[i] <- auc_value_train
  }
  avg_aucs[model_name] <- mean(auc_values_train)
}

# Convert the average AUCs to a matrix
avg_aucs_matrix <- matrix(avg_aucs, nrow = 1, dimnames = list("Average AUC", names(avg_aucs)))
print(avg_aucs_matrix)

###############################################################################
############################### Overdispersion 1 ############################## 
###############################################################################

# Function to calculate overdispersion for each model
calculate_overdispersion <- function(model) {
  residual_deviance <- deviance(model)
  residual_degrees_of_freedom <- df.residual(model)
  dispersion_ratio <- residual_deviance / residual_degrees_of_freedom
  return(dispersion_ratio)
}

# Extract component models from the model averaging object
component_models <- get.models(avg_list_feed_30m[[1]], subset = TRUE)

# Calculate overdispersion for each component model
dispersion_ratios <- sapply(component_models, calculate_overdispersion)

# Print dispersion ratios for each model
print(dispersion_ratios)

# Calculate mean dispersion ratio
mean_dispersion_ratio <- mean(dispersion_ratios, na.rm = TRUE)

# Print mean dispersion ratio
print(paste("Mean Dispersion Ratio: ", mean_dispersion_ratio))

# Check if the averaged model is overdispersed
if (mean_dispersion_ratio > 2) {
  print("The model is overdispersed")
} else if (mean_dispersion_ratio > 1) {
  print("The model may be overdispersed")
} else {
  print("The model is not overdispersed")
}

###############################################################################
############################## Overdispersion 100 ############################# 
###############################################################################

load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/glm_avg_list/avg_list_nest_scale.RData")

# Function to calculate overdispersion for all models in the list
calculate_overdispersion_for_all <- function(avg_list) {
  overdispersion_list <- list()
  
  for (i in 1:length(avg_list)) {
    avg_model <- avg_list[[i]]
    component_models <- get.models(avg_model, subset = TRUE)
    dispersion_ratios <- sapply(component_models, calculate_overdispersion)
    mean_dispersion_ratio <- mean(dispersion_ratios, na.rm = TRUE)
    overdispersion_list[[i]] <- mean_dispersion_ratio
  }
  
  return(overdispersion_list)
}

# Calculate overdispersion for all models in avg_list_feed_30m
overdispersion_results <- calculate_overdispersion_for_all(avg_list_feed_30m)

# Print overdispersion results for each model
print(overdispersion_results)

# Calculate summary statistics for overdispersion results
overdispersion_summary <- data.frame(
  Metric = "Dispersion Ratio",
  Min = min(unlist(overdispersion_results), na.rm = TRUE),
  Max = max(unlist(overdispersion_results), na.rm = TRUE),
  Mean = mean(unlist(overdispersion_results), na.rm = TRUE),
  SD = sd(unlist(overdispersion_results), na.rm = TRUE)
)

# Print summary statistics
print(overdispersion_summary)

# Check if any of the averaged models are overdispersed
overdispersed_models <- sum(unlist(overdispersion_results) > 2)
potentially_overdispersed_models <- sum(unlist(overdispersion_results) > 1)

print(paste("Number of overdispersed models: ", overdispersed_models))
print(paste("Number of potentially overdispersed models: ", potentially_overdispersed_models))

###############################################################################
############################### Models Overlap ################################
###############################################################################

load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_scale.RData")

load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/glm/V2/predictions_df_glm_nest_scale.RData")

library(ggplot2)
library(raster)
library(sp)
library(sf)
library(fields)

# Load the MCP boundary
elev_30m_MCP <- crop(elev_30m, extent(MCP_buffered))
load("/Users/lucasdegreve/Data/RData/centroids_plot.RData")
centroids_sf <- st_as_sf(centroids_df, coords = c("x", "y"), crs = st_crs(Nimba_WHS_buffered))
MCP_buffered <- st_transform(MCP_buffered, st_crs(centroids_sf))
inside_MCP <- st_within(centroids_sf, MCP_buffered, sparse = FALSE)
centroids_MCP <- centroids_plot
centroids_MCP[!inside_MCP, c("x", "y")] <- NA
Nimba_WHS_within_MCP <- st_intersection(Nimba_WHS, MCP_buffered)
valid_indices <- !is.na(centroids_MCP$x) & !is.na(centroids_MCP$y)
filtered_centroids_MCP <- centroids_MCP[valid_indices, ]

# Define plot dimensions
plot_extent <- extent(elev_30m_MCP)
plot_width <- plot_extent@xmax - plot_extent@xmin
plot_height <- plot_extent@ymax - plot_extent@ymin

# Function to generate and save plots
generate_plots <- function(predictions1, predictions2, predictions3, model_name) {
  prediction_means1 <- rowMeans(predictions1)
  prediction_means2 <- rowMeans(predictions2)
  prediction_means3 <- rowMeans(predictions3)
  
  prediction_means <- prediction_means1*prediction_means2
  plot_file1 <- paste0("plot_", model_name, "_feed.png")
  
  png(plot_file1, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  filtered_prediction_means <- prediction_means[valid_indices]
  min_means <- min(filtered_prediction_means)
  max_means <- max(filtered_prediction_means)
  colourScale = colorRampPalette(brewer.pal(11, "OrRd"))(121)[11:121]
  color_indices <- (((prediction_means - min_means) / (max_means - min_means)) * 100) + 1
  color_indices <- pmax(1, pmin(100, color_indices))
  cols <- colourScale[color_indices]
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", frame.plot = FALSE)
  image.plot(z = filtered_prediction_means, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
  
  prediction_means <- prediction_means3*prediction_means2
  plot_file2 <- paste0("plot_", model_name, "_nest.png")
  
  png(plot_file2, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  filtered_prediction_means <- prediction_means[valid_indices]
  min_means <- min(filtered_prediction_means)
  max_means <- max(filtered_prediction_means)
  colourScale = colorRampPalette(brewer.pal(11, "OrRd"))(121)[11:121]
  color_indices <- (((prediction_means - min_means) / (max_means - min_means)) * 100) + 1
  color_indices <- pmax(1, pmin(100, color_indices))
  cols <- colourScale[color_indices]
  plot(centroids_MCP$x, centroids_MCP$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", frame.plot = FALSE)
  image.plot(z = filtered_prediction_means, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_within_MCP, border = "black", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "white", col = "transparent", add = TRUE)
  plot(MCP_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered survey area"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
}

# Apply the function to each model combination
generate_plots(prediction_df_brt_feed_scale, prediction_df_brt_hunt_scale, prediction_df_brt_nest_scale, "brt_scale")
generate_plots(prediction_df_brt_feed_30m, prediction_df_brt_hunt_30m, prediction_df_brt_nest_30m, "brt_30m")
generate_plots(predictions_df_glm_feed_30m, predictions_df_glm_hunt_30m, predictions_df_glm_nest_30m, "glm_30m")
generate_plots(predictions_df_glm_feed_scale, predictions_df_glm_hunt_scale, predictions_df_glm_nest_scale, "glm_scale")
