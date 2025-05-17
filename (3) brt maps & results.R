###############################################################################
################################### loading ###################################
###############################################################################

library(sf)
library(readr)
library(raster)
library(fields)
library(RColorBrewer)
library(gbm)
library(openxlsx)
library(sp)

# Load the tracks and the reserve boundary
tracks <- read_sf("Data/tracks/tracks.shp")
MCP <- read_sf("Data/MCP/MCP.shp")
Nimba_WHS <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS.shp")
Nimba_WHS_buffered <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS_buffered.shp")

# Load the presence points
nesting <- read_delim("Data/presence/nesting.csv", delim = ",")
hunting <- read_delim("Data/presence/hunting.csv", delim = ",")
feeding_travelling <- read_delim("Data/presence/feeding-travelling.csv", delim = ",")
feed_coords = st_as_sf(feeding_travelling, coords=c("X","Y"))
nest_coords = st_as_sf(nesting, coords=c("X","Y"))
hunt_coords = st_as_sf(hunting, coords=c("X","Y"))

# Load elevation raster and centroids
elev_30m <- raster("GEE/30m/elev_30m.tif")
elev_30m <- crop(elev_30m, extent(Nimba_WHS_buffered))
centroids <- rasterToPoints(elev_30m, spatial = TRUE)
centroids_df <- data.frame(centroids@coords)
load("/Users/lucasdegreve/Data/RData/centroids_plot.RData")

load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_scale.RData")

###############################################################################
################################## centroids ##################################
###############################################################################

elev_tracks <- extract(elev_30m, as(tracks, "Spatial"))
elev_tracks <- unlist(elev_tracks)
max_elev_tracks <- max(elev_tracks, na.rm = TRUE)
elev_values <- extract(elev_30m, centroids_df)
centroids_sf <- st_as_sf(centroids_df, coords = c("x", "y"), crs = st_crs(Nimba_WHS_buffered))
within_reserve <- st_within(centroids_sf, Nimba_WHS_buffered, sparse = FALSE)
centroids_plot <- centroids_df # Centroids for plotting
centroids_plot[elev_tracks > max_elev_tracks | !within_reserve, ] <- NA
save(centroids_plot, file = "centroids_plot.RData")

###############################################################################
################################## BRT maps ###################################
###############################################################################

plot_prediction_means <- function(prediction_df, centroids_plot, plot_width, plot_height, Nimba_WHS, Nimba_WHS_buffered, output_file) {
  png(output_file, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  
  prediction_means <- rowMeans(prediction_df)
  min_means <- min(prediction_means)
  max_means <- max(prediction_means)
  colourScale = colorRampPalette(brewer.pal(11, "YlGn"))(121)[11:121]
  cols <- colourScale[(((prediction_means - min_means) / (max_means - min_means)) * 100) + 1]
  
  plot(centroids_plot$x, centroids_plot$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", main = "Prediction", frame.plot = FALSE)
  image.plot(z = prediction_means, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  plot(Nimba_WHS, border = "black", col = "transparent", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered reserve"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
}

plot_prediction_sd <- function(prediction_df, centroids_plot, plot_width, plot_height, Nimba_WHS, Nimba_WHS_buffered, output_file) {
  png(output_file, width = plot_width, height = plot_height, units = "px", res = 1900)
  par(mar = c(9, 4, 1, 1))  # c(bottom, left, top, right)
  
  prediction_sd <- apply(prediction_df, 1, sd)
  min_sd <- min(prediction_sd)
  max_sd <- max(prediction_sd)
  colourScale <- colorRampPalette(brewer.pal(11, "OrRd"))(121)[11:121]
  cols <- colourScale[(((prediction_sd - min_sd) / (max_sd - min_sd)) * 100) + 1]
  
  plot(centroids_plot$x, centroids_plot$y, pch = ".", col = cols, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", main = "Standard deviation", frame.plot = FALSE) # Plot the map
  image.plot(z = prediction_sd, col = colourScale, legend.only = TRUE, horizontal = TRUE, axis.args = list(at = NULL, labels = NULL))
  plot(Nimba_WHS_buffered, border = "black", col = "transparent", lty = "dashed", add = TRUE)
  plot(Nimba_WHS, border = "black", col = "transparent", add = TRUE)
  legend("bottomright", legend = c("Mount Nimba Strict Nature Reserve", "2-kilometer buffered reserve"), col = c("black", "black"), lty = c(1, 2), cex = 0.8, bty = "n")
  dev.off()  # Close PNG device
}

# List of prediction datasets
prediction_dfs <- list(
  prediction_df_brt_feed_30m = prediction_df_brt_feed_30m,
  prediction_df_brt_feed_scale = prediction_df_brt_feed_scale,
  prediction_df_brt_hunt_30m = prediction_df_brt_hunt_30m,
  prediction_df_brt_hunt_scale = prediction_df_brt_hunt_scale,
  prediction_df_brt_nest_30m = prediction_df_brt_nest_30m,
  prediction_df_brt_nest_scale = prediction_df_brt_nest_scale
)

# Define plot dimensions
plot_extent <- extent(elev_30m)
plot_width <- plot_extent@xmax - plot_extent@xmin
plot_height <- plot_extent@ymax - plot_extent@ymin

# Loop through each prediction dataframe and create plots
for (name in names(prediction_dfs)) {
  prediction_df <- prediction_dfs[[name]]
  # Create prediction means plot
  plot_prediction_means(prediction_df, centroids_plot, plot_width, plot_height, 
    Nimba_WHS, Nimba_WHS_buffered, paste0(name, "_prediction_means.png")
  )
  # Create prediction standard deviation plot
  plot_prediction_sd(prediction_df, centroids_plot, plot_width, plot_height, 
    Nimba_WHS, Nimba_WHS_buffered, paste0(name, "_prediction_sd.png")
  )
}

###############################################################################
##################################### AUC #####################################
###############################################################################

load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_scale_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_30m_list_brt_results.RData")
load("/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_scale_list_brt_results.RData")

extract_auc_summary <- function(brt_list) {
  auc_values <- sapply(brt_list, function(brt_model) brt_model$cv.statistics$discrimination.mean)
  auc_summary <- c(mean = mean(auc_values), min = min(auc_values), max = max(auc_values), sd = sd(auc_values))
  return(auc_summary)
}

auc_brt_feed_30m <- extract_auc_summary(df_brt_feed_30m_list_brt_results)
auc_brt_feed_scale <- extract_auc_summary(df_brt_feed_scale_list_brt_results)
auc_brt_hunt_30m <- extract_auc_summary(df_brt_hunt_30m_list_brt_results)
auc_brt_hunt_scale <- extract_auc_summary(df_brt_hunt_scale_list_brt_results)
auc_brt_nest_30m <- extract_auc_summary(df_brt_nest_30m_list_brt_results)
auc_brt_nest_scale <- extract_auc_summary(df_brt_nest_scale_list_brt_results)

auc_table <- data.frame(
  Model = c("feed_30m", "feed_scale", "hunt_30m", "hunt_scale", "nest_30m", "nest_scale"),
  Mean = c(auc_brt_feed_30m[1], auc_brt_feed_scale[1], auc_brt_hunt_30m[1], auc_brt_hunt_scale[1], auc_brt_nest_30m[1], auc_brt_nest_scale[1]),
  Min  = c(auc_brt_feed_30m[2], auc_brt_feed_scale[2], auc_brt_hunt_30m[2], auc_brt_hunt_scale[2], auc_brt_nest_30m[2], auc_brt_nest_scale[2]),
  Max  = c(auc_brt_feed_30m[3], auc_brt_feed_scale[3], auc_brt_hunt_30m[3], auc_brt_hunt_scale[3], auc_brt_nest_30m[3], auc_brt_nest_scale[3]),
  SD   = c(auc_brt_feed_30m[4], auc_brt_feed_scale[4], auc_brt_hunt_30m[4], auc_brt_hunt_scale[4], auc_brt_nest_30m[4], auc_brt_nest_scale[4])
)

write.xlsx(auc_table, file = "brt_auc.xlsx", sheetName = "BRT AUC", colNames = TRUE, rowNames = FALSE)

###############################################################################
############################# Relative Influences #############################
###############################################################################

library(gbm)

load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_feed_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_feed_scale_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_hunt_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_hunt_scale_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_nest_30m_list.RData")
load("/Users/lucasdegreve/Data/RData/df_brt/df_brt_nest_scale_list.RData")

# Define the covariates for each model
covariates <- list(
  brt_feed_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_feed_scale = c("aspect_1500k", "cti_1500k", "elev_30m", "evi_100k", "hli_2000k", "slope_250k", "tpi_2000k", "treecover_2000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_hunt_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_hunt_scale = c("aspect_2000k", "cti_1000k", "elev_30m", "evi_1500k", "hli_2000k", "slope_100k", "tpi_1000k", "treecover_500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_nest_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_nest_scale = c("aspect_1500k", "cti_500k", "elev_100k", "evi_1000k", "hli_1500k", "slope_500k", "tpi_2000k", "treecover_2000k", "wetness_2000k", "distance_camp", "distance_river", "distance_road", "distance_village")
)

# Store the list of BRT results in a named list
brt_results <- list(
  df_brt_feed_30m_list_brt_results = df_brt_feed_30m_list_brt_results,
  df_brt_feed_scale_list_brt_results = df_brt_feed_scale_list_brt_results,
  df_brt_hunt_30m_list_brt_results = df_brt_hunt_30m_list_brt_results,
  df_brt_hunt_scale_list_brt_results = df_brt_hunt_scale_list_brt_results,
  df_brt_nest_30m_list_brt_results = df_brt_nest_30m_list_brt_results,
  df_brt_nest_scale_list_brt_results = df_brt_nest_scale_list_brt_results
)

# Function to calculate relative influences
calculate_relative_influences <- function(brt_list, covariates) {
  relativeInfluences <- matrix(0, nrow=length(covariates), ncol=1)
  for (j in 1:length(brt_list)) {
    model_summary <- summary(brt_list[[j]])
    for (k in 1:length(covariates)) {
      covariate <- covariates[k]
      if (covariate %in% rownames(model_summary)) {
        relativeInfluences[k] <- relativeInfluences[k] + model_summary[covariate, "rel.inf"]
      }
    }
  }
  relativeInfluences <- relativeInfluences / length(brt_list)
  row.names(relativeInfluences) <- covariates
  return(relativeInfluences)
}

# Function to create response curves
create_response_curves <- function(brt_list, covariates, data, relativeInfluences, output_filename) {
  envVariableValues <- matrix(nrow=3, ncol=length(covariates))
  row.names(envVariableValues) <- c("median", "minV", "maxV")
  colnames(envVariableValues) <- covariates
  for (j in 1:length(covariates)) {
    minV <- min(data[, covariates[j]], na.rm=TRUE)
    maxV <- max(data[, covariates[j]], na.rm=TRUE)
    medianV <- median(data[, covariates[j]], na.rm=TRUE)
    envVariableValues[, j] <- c(medianV, minV, maxV)
  }
  envVariableValues_list <- envVariableValues
  pdf(output_filename, width = 8.75, height = 5.30)
  par(mfrow=c(3,5), mar=c(2,1,0.5,1), lwd=1, col="black")
  for (i in 1:length(covariates)) {
    predictions_list <- list()
    valuesInterval <- (envVariableValues_list["maxV",i] - envVariableValues_list["minV",i]) / 100
    df <- data.frame(matrix(nrow=length(seq(envVariableValues_list["minV",i], envVariableValues_list["maxV",i], valuesInterval)), ncol=length(covariates)))
    colnames(df) <- covariates
    for (k in 1:length(covariates)) {
      valuesInterval <- (envVariableValues_list["maxV",k] - envVariableValues_list["minV",k]) / 100
      if (i == k) {
        df[, covariates[k]] <- seq(envVariableValues_list["minV", k], envVariableValues_list["maxV", k], valuesInterval)
      } else {
        df[, covariates[k]] <- rep(envVariableValues_list["median", k], nrow(df))
      }
    }
    for (j in 1:length(brt_list)) {
      n.trees <- brt_list[[j]]$gbm.call$best.trees
      prediction <- predict.gbm(brt_list[[j]], newdata=df, n.trees, type="response")
      predictions_list[[j]] <- prediction
    }
    plot(df[, covariates[i]], predictions_list[[1]], col="red", ann=FALSE, axes=FALSE, lwd=0.2, type="l")
    for (l in 2:length(predictions_list)) {
      lines(df[, covariates[i]], predictions_list[[l]], col="red", lwd=0.2)
    }
    axis(side=1, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.07,0))
    axis(side=2, lwd.tick=0.2, cex.axis=0.7, lwd=0, tck=-0.030, col.axis="gray30", mgp=c(0,0.2,0))
    title(ylab="predicted values", cex.lab=0.9, mgp=c(1.3,0,0), col.lab="gray30")
    title(xlab=paste0(covariates[i], " (", round(relativeInfluences[i, 1], 1), "%)"), cex.lab=0.9, mgp=c(0.9,0,0), col.lab="gray30")
    box(lwd=0.2, col="gray30")
  }
  dev.off()
}

relative_influences_list <- list()
# Iterate over each model name
for (model_name in names(covariates)) {
  covariate_set <- covariates[[model_name]]
  brt_list <- brt_results[[paste0("df_", model_name, "_list_brt_results")]]
  data_list <- get(paste0("df_", model_name, "_list"))
  relative_influences <- calculate_relative_influences(brt_list, covariate_set)
  relative_influences_list[[model_name]] <- relative_influences
  data_curv <- data_list[[1]][which(data_list[[1]][, "presence"] == 1),]
  output_filename <- paste0(model_name, "_response_curves.pdf")
  create_response_curves(brt_list, covariate_set, data_curv, relative_influences, output_filename)
}
# Convert relative_influences_list to a data frame
relative_influences_table <- as.data.frame(relative_influences_list)

library(openxlsx)
write.xlsx(relative_influences_table, file = "relative_influences.xlsx", sheetName = "Relative Influences", colNames = TRUE, rowNames = TRUE)

###############################################################################
################################## Box Plot ################################### 
###############################################################################

# Initialize matrices and data frames for storing results
relativeInfluences_all <- list()
relativeInfluences_all_df <- data.frame()
relativeInfluences_all_transposed <- list()

# Loop through each model
for (model_name in names(brt_results)) {
  cat("Processing model:", model_name, "\n")
  # Get covariates and BRT results for the current model
  covariate_set <- covariates[[gsub("df_", "", gsub("_list_brt_results", "", model_name))]]
  brt_list <- brt_results[[model_name]]
  data_list <- get(gsub("_list_brt_results", "_list", model_name))
  # Calculate relative influences for the current model
  relativeInfluences_all[[model_name]] <- matrix(0, nrow=length(covariate_set), ncol=length(brt_list))
  row.names(relativeInfluences_all[[model_name]]) <- covariate_set
  for (i in 1:length(brt_list)) {
    for (k in 1:length(covariate_set)) {
      relativeInfluences_all[[model_name]][k, i] <- summary(brt_list[[i]])[gsub("-", "\\.", covariate_set)[k], "rel.inf"]
    }
  }
  # Compute means and confidence intervals
  RI_Means <- rowMeans(relativeInfluences_all[[model_name]])
  ci_lower <- c()
  ci_upper <- c()
  for (i in 1:length(covariate_set)) {
    # Calculate 95% confidence interval
    ci <- t.test(relativeInfluences_all[[model_name]][i, ])$conf.int
    # Append results to lists
    ci_lower <- c(ci_lower, ci[1])
    ci_upper <- c(ci_upper, ci[2])
  }
  # Create a data frame for the current model
  model_df <- data.frame(
    Model = rep(gsub("df_", "", gsub("_list_brt_results", "", model_name)), length(covariate_set)),
    Covariate = rep(covariate_set, each = length(brt_list)),
    RI_Means = RI_Means,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  )
  # Append to the main data frame
  relativeInfluences_all_df <- rbind(relativeInfluences_all_df, model_df)
  # Transpose relative influences matrix for boxplot
  relativeInfluences_all_transposed[[model_name]] <- t(relativeInfluences_all[[model_name]])

  # Boxplot for the current model
  pdf(paste0(gsub("df_", "", gsub("_list_brt_results", "", model_name)), "_boxplot.pdf"), width = 20, height = 10)
  par(mfrow = c(1, 1), mar = c(3, 3, 3, 3), lwd = 1, col = "black")
  boxplot(relativeInfluences_all_transposed[[model_name]], main = model_name, xlab = "Covariates", ylab = "Relative Influence")
  dev.off()
}
write.xlsx(relativeInfluences_all_df, file = "relative_influences_all_models.xlsx", sheetName = "Relative Influences", colNames = TRUE, rowNames = FALSE)

###############################################################################
################################## Sorensen ################################### 
###############################################################################

load("/Users/lucasdegreve/Data/RData/pa_list/pa_feed_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_hunt_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_nest_list.RData")

Nimba_WHS_buffered <- read_sf("QGIS/Layers for Lucas/Mount Nimba Strict Nature Reserve_boundary/Nimba_WHS_buffered.shp")
elev_30m <- raster("GEE/30m/elev_30m.tif")
elev_30m <- crop(elev_30m, extent(Nimba_WHS_buffered))
centroids <- rasterToPoints(elev_30m, spatial = TRUE)
centroids_df <- data.frame(centroids@coords)

# Load prediction data frames
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_feed_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_hunt_scale.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_30m.RData")
load("/Users/lucasdegreve/Data/RData/prediction_df/prediction_df_brt_nest_scale.RData")

models <- list(
  "feed_30m" = list(prediction_df = "prediction_df_brt_feed_30m", pa_list = "pa_feed_list"),
  "feed_scale" = list(prediction_df = "prediction_df_brt_feed_scale", pa_list = "pa_feed_list"),
  "hunt_30m" = list(prediction_df = "prediction_df_brt_hunt_30m", pa_list = "pa_hunt_list"),
  "hunt_scale" = list(prediction_df = "prediction_df_brt_hunt_scale", pa_list = "pa_hunt_list"),
  "nest_30m" = list(prediction_df = "prediction_df_brt_nest_30m", pa_list = "pa_nest_list"),
  "nest_scale" = list(prediction_df = "prediction_df_brt_nest_scale", pa_list = "pa_nest_list")
)

# Function to process each model
process_model <- function(model_name, prediction_df, pa_list) {
  prediction_coords <- cbind(centroids_df, get(prediction_df))
  data_sorensen_list <- list()
  
  for (i in 1:100) {
    pa <- get(pa_list)[[i]][[1]]  # Adjusting to the list structure
    coord <- get(pa_list)[[i]][[2]]
    pa_sp <- SpatialPointsDataFrame(coords = coord, data = pa)
    coordinates_pa <- coordinates(pa_sp)
    coordinates_pred <- as.matrix(prediction_coords[, 1:2])
    match_indices <- match(paste(coordinates_pa[, 1], coordinates_pa[, 2]), paste(coordinates_pred[, 1], coordinates_pred[, 2]))
    predictions <- prediction_coords[match_indices, i + 2]  # +2 to skip 'x' and 'y' columns in prediction_coords
    data_sorensen <- data.frame(presence = pa$presence, prediction = predictions)
    data_sorensen_list[[i]] <- data_sorensen
  }
  data_comp_list <- list()
  
  for (k in 1:100) {
    sorensen_df <- data_sorensen_list[[k]]
    thres <- seq(0.01, 1, by = 0.01)
    data_comp <- as.data.frame(matrix(0, nrow = length(thres), ncol = 6))
    colnames(data_comp) <- c("TP", "FP", "FN", "TN", "SI", "threshold")
    data_comp[,"threshold"] <- thres
    
    for (i in seq_along(thres)) {
      for (j in 1:nrow(sorensen_df)) {
        # True positive
        if (sorensen_df[j, 1] == 1 && sorensen_df[j, 2] >= thres[i]) {
          data_comp[i, "TP"] <- data_comp[i, "TP"] + 1
        }
        # False negative
        if (sorensen_df[j, 1] == 1 && sorensen_df[j, 2] < thres[i]) {
          data_comp[i, "FN"] <- data_comp[i, "FN"] + 1
        }
        # False positive
        if (sorensen_df[j, 1] == 0 && sorensen_df[j, 2] >= thres[i]) {
          data_comp[i, "FP"] <- data_comp[i, "FP"] + 1
        }
        # True negative
        if (sorensen_df[j, 1] == 0 && sorensen_df[j, 2] < thres[i]) {
          data_comp[i, "TN"] <- data_comp[i, "TN"] + 1
        }
      }
    }
    x <- sum(sorensen_df$presence == 1) / sum(sorensen_df$presence == 0) * (1 - sum(sorensen_df$presence == 1) / nrow(sorensen_df)) / (sum(sorensen_df$presence == 1) / nrow(sorensen_df))
    data_comp[,"SI"] <- (2 * data_comp[,"TP"]) / ((2 * data_comp[,"TP"]) + (x * data_comp[,"FP"]) + data_comp[,"FN"])
    data_comp_list[[k]] <- data_comp
  }
  max_si_df <- data.frame(matrix(ncol = 2, nrow = 100))
  colnames(max_si_df) <- c("Set", "Max_SI")
  
  for (i in 1:100) {
    data_comp <- data_comp_list[[i]]
    max_si <- max(data_comp$SI, na.rm = TRUE)
    max_si_df[i, ] <- c(i, max_si)
  }
  si_means <- mean(max_si_df$Max_SI, na.rm = TRUE)
  si_sd <- sd(max_si_df$Max_SI, na.rm = TRUE)
  print(paste("Mean Sorensen Index for", model_name, ":", si_means))
  print(paste("SD of Sorensen Index for", model_name, ":", si_sd))
  write.xlsx(data_comp, file = paste0("data_comp_", model_name, ".xlsx"), sheetName = paste("Sorensen", model_name), colNames = TRUE, rowNames = TRUE)
}

# Loop through each model and process it
for (model_name in names(models)) {
  prediction_df <- models[[model_name]]$prediction_df
  pa_list <- models[[model_name]]$pa_list
  process_model(model_name, prediction_df, pa_list)
}

###############################################################################
################################### Env Var ################################### 
###############################################################################

par(mfrow=c(4,4), mar=c(1,1,1,1), lwd=1)
for (variable_name in brt_feed_30m) {
  raster_name <- ifelse(variable_name %in% c("camp", "river", "road", "village"),
                        paste0("distance_", variable_name),
                        paste0(variable_name, "_30m"))
  plot(get(raster_name), main = "")
  title(main = variable_name, line = +1)  # Adjust the line parameter to move the title up
}

tracks <- read_sf("QGIS/Nimba GIS Data/clean tracks/tracks_off.shp")
plot(elev_30m)
plot(tracks, col = "black", add = TRUE)


png("plot.png", width = 24900, height = 25300, units = "px", res = 1900)
colourScale = colorRampPalette(brewer.pal(11, "YlOrBr"))(121)[11:121]
plot(elev_30m, col = colourScale, cex = 1, lwd = 1, xlab = "Longitude", ylab = "Latitude", frame.plot = FALSE)
plot(Nimba_WHS, border = "white", col = "transparent", add = TRUE)
plot(tracks, col = "black", add = TRUE)
legend("bottomright", legend = c("Tracks surveyed", "Mount Nimba Strict Nature Reserve"), col = c("black", "black"), lty = c(1, 2), cex = 0.8)
dev.off()  # Close PNG device

tracks_raster <- raster("Data/tracks_raster.tif")
tracks_raster <- crop(tracks_raster, extent(MCP))
plot(tracks_raster, colNA = "black")
