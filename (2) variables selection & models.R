###############################################################################
################################### loading ###################################
###############################################################################

# Set directory
setwd("/Users/lucasdegreve/")

# Load necessary libraries
library(readr)
library(sf)
library(raster)
library(pROC)
library(sp)
library(dismo)
library(blockCV)
library(Hmisc)
library(corrplot)
library(MuMIn)
library(gbm)

# Load the presence points
nesting <- read_delim("Data/presence/nesting.csv", delim = ",")
feeding_travelling <- read_delim("Data/presence/feeding-travelling.csv", delim = ",")
hunting <- read_delim("Data/presence/hunting.csv", delim = ",")

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

# Load the presence-absence lists
load("/Users/lucasdegreve/Data/RData/pa_list/pa_feed_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_hunt_list.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_nest_list.RData")

###############################################################################
################################## dataframes #################################
###############################################################################


# Define the variables and scales
distance_variables <- c("distance_camp", "distance_river", "distance_road", "distance_village")
scales <- c("30m", "100k", "250k", "500k", "1000k", "1500k", "2000k")
variables <- c("aspect", "brightness", "cti", "elev", "evi", "greenness", "hli", "roughness", "slope", "tpi", "treecover", "wetness")

# Function to create data frames for a given list and type (feed, hunt, nest)
create_data_frames <- function(pa_list, type) {
  df_list <- list()
  for (i in 1:100) {
    pa_data <- pa_list[[i]][[paste0("pa_", type)]]
    coord_data <- pa_list[[i]][[paste0("coord_", type)]]
    scale_list <- list()
    variable_names <- c("x", "y", "presence")
    for (variable in variables) {
      for (scale in scales) {
        raster_data <- extract(get(paste(variable, scale, sep = "_")), coord_data)
        scale_list[[length(scale_list) + 1]] <- raster_data
        variable_names <- c(variable_names, paste(variable, scale, sep = "_"))
      }
    }
    for (distance_var in distance_variables) {
      raster_data <- extract(get(distance_var), coord_data)
      scale_list[[length(scale_list) + 1]] <- raster_data
      variable_names <- c(variable_names, distance_var)
    }
    df <- cbind(pa_data, do.call(cbind, scale_list))
    colnames(df) <- variable_names
    df_list[[i]] <- df
  }
  return(df_list)
}

# Create data frames for each list
df_feed_list <- create_data_frames(pa_feed_list, "feed")
df_hunt_list <- create_data_frames(pa_hunt_list, "hunt")
df_nest_list <- create_data_frames(pa_nest_list, "nest")

save(df_feed_list, file = "df_feed_list.RData")
save(df_hunt_list, file = "df_hunt_list.RData")
save(df_nest_list, file = "df_nest_list.RData")

# Function to apply create_data_frames to new 85 lists
apply_create_data_frames <- function(list85, type) {
  df_85_list <- create_data_frames(list85, paste0(type, "85"))
  return(df_85_list)
}

# Create data frames for each of the new 85 lists
df_feed_list85 <- apply_create_data_frames(pa_feed_list85, "feed")
df_hunt_list85 <- apply_create_data_frames(pa_hunt_list85, "hunt")
df_nest_list85 <- apply_create_data_frames(pa_nest_list85, "nest")

###############################################################################
############################### Uni-variate BRT ###############################
###############################################################################

# Load necessary libraries
library(dismo)
library(blockCV)

# Function to perform BRT and calculate AUC for a given dataset list
calculate_auc <- function(df_list, n.folds = 5, theRanges = c(1300, 1300)) {
  auc_matrix <- matrix(NA, nrow = 100, ncol = 88)
  for (i in 1:100) {
    dataset_i <- df_list[[i]]
    spdf <- SpatialPointsDataFrame(coords = dataset_i[, c("x", "y")], data = dataset_i, proj4string = crs(elev_30m))
    # Generate cross-validation folds
    folds_with_similar_sizes <- FALSE
    while (!folds_with_similar_sizes) {
      myblocks <- cv_spatial(spdf, column ="presence", k=n.folds, size=theRanges, selection="random", progress = FALSE)
      fold.vector3 <- myblocks$folds_ids
      fold.vector_presences <- fold.vector3[which(dataset_i[, "presence"] == 1)]
      counts <- hist(fold.vector_presences, plot = FALSE)$counts
      props <- counts[which(counts > 0)] / sum(counts)
      if (min(props) > 0.05) folds_with_similar_sizes <- TRUE
    }
    # Loop over each predictor variable (starting from column 4)
    for (j in 1:88) {
      gbm.x <- j + 3
      # Fit BRT model
      brt_model <- gbm.step(dataset_i, gbm.x, gbm.y = 3, offset = NULL, fold.vector = fold.vector3,
                            tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.80,
                            site.weights = rep(1, dim(dataset_i)[1]), var.monotone = rep(0, length(gbm.x)), n.folds = 5,
                            prev.stratify = TRUE, family = "bernoulli", n.trees = 100, step.size = 10, max.trees = 10000,
                            tolerance.method = "auto", tolerance = 0.001, plot.main = TRUE, plot.folds = FALSE, verbose = TRUE,
                            silent = FALSE, keep.fold.models = FALSE, keep.fold.vector = FALSE, keep.fold.fit = FALSE)
      if (!is.null(brt_model$cv.statistics$discrimination.mean)) {
        auc_matrix[i, j] <- brt_model$cv.statistics$discrimination.mean
      }
    }
  }
  return(auc_matrix)
}

# Calculate AUC for each list
auc_feed <- calculate_auc(df_feed_list)
auc_hunt <- calculate_auc(df_hunt_list)
auc_nest <- calculate_auc(df_nest_list)

# Rename columns of the AUC matrices
colnames(auc_feed) <- colnames(df_feed_list[[1]])[4:91] # Extract column names from the df_feed_list (columns 4 to 91)
colnames(auc_hunt) <- colnames(df_feed_list[[1]])[4:91]
colnames(auc_nest) <- colnames(df_feed_list[[1]])[4:91]

# Save each AUC matrix as a CSV file
write.csv(auc_feed, "auc_brt_feed.csv", row.names = FALSE)
write.csv(auc_hunt, "auc_brt_hunt.csv", row.names = FALSE)
write.csv(auc_nest, "auc_brt_nest.csv", row.names = FALSE)

auc_brt_feed <- read.csv("~/Data/auc & aic/auc_brt_feed.csv")
auc_brt_hunt <- read.csv("~/Data/auc & aic/auc_brt_hunt.csv")
auc_brt_nest <- read.csv("~/Data/auc & aic/auc_brt_nest.csv")

###############################################################################
############################### Uni-variate GLM ###############################
###############################################################################

load("/Users/lucasdegreve/Data/RData/df_list/df_feed_list.RData")
load("/Users/lucasdegreve/Data/RData/df_list/df_hunt_list.RData")
load("/Users/lucasdegreve/Data/RData/df_list/df_nest_list.RData")

library(pROC)

# Function to calculate AUC and AIC for a given list
calculate_auc_aic <- function(df_list) {
  auc_glm <- matrix(NA, nrow = 100, ncol = 88)
  aic_glm <- matrix(NA, nrow = 100, ncol = 88)
  for (i in 1:100) { # Loop over each dataset in df_list
    dataset_glm <- df_list[[i]]
    for (j in 1:88) { # Loop over each predictor variable
      glm_model <- glm(presence ~ dataset_glm[, j+3], data = dataset_glm, family = binomial(link = "logit"))   # Fit univariate GLM
      pred_probs <- predict(glm_model, type = "response")   # Predict probabilities
      auc_glm[i, j] <- roc(dataset_glm$presence, pred_probs)$auc
      aic_glm[i, j] <- AIC(glm_model)
    }
  }
  return(list(auc_glm = auc_glm, aic_glm = aic_glm))
}

# Apply the function to each list
result_feed <- calculate_auc_aic(df_feed_list)
result_hunt <- calculate_auc_aic(df_hunt_list)
result_nest <- calculate_auc_aic(df_nest_list)

# Extract the results
auc_glm_feed <- result_feed$auc_glm
aic_glm_feed <- result_feed$aic_glm
auc_glm_hunt <- result_hunt$auc_glm
aic_glm_hunt <- result_hunt$aic_glm
auc_glm_nest <- result_nest$auc_glm
aic_glm_nest <- result_nest$aic_glm

colnames(auc_glm_feed) <- colnames(df_feed_list[[1]])[4:91] # Extract column names from the df_feed_list (columns 4 to 91)
colnames(aic_glm_feed) <- colnames(df_feed_list[[1]])[4:91]
colnames(auc_glm_hunt) <- colnames(df_feed_list[[1]])[4:91]
colnames(aic_glm_hunt) <- colnames(df_feed_list[[1]])[4:91]
colnames(auc_glm_nest) <- colnames(df_feed_list[[1]])[4:91]
colnames(aic_glm_nest) <- colnames(df_feed_list[[1]])[4:91]

# Save each AUC matrix as a CSV file
write.csv(auc_glm_feed, "auc_glm_feed.csv", row.names = FALSE)
write.csv(aic_glm_feed, "aic_glm_feed.csv", row.names = FALSE)
write.csv(auc_glm_hunt, "auc_glm_hunt.csv", row.names = FALSE)
write.csv(aic_glm_hunt, "aic_glm_hunt.csv", row.names = FALSE)
write.csv(auc_glm_nest, "auc_glm_nest.csv", row.names = FALSE)
write.csv(aic_glm_nest, "aic_glm_nest.csv", row.names = FALSE)

auc_brt_feed <- read.csv("~/Data/auc & aic/auc_brt_feed.csv")
auc_brt_hunt <- read.csv("~/Data/auc & aic/auc_brt_hunt.csv")
auc_brt_nest <- read.csv("~/Data/auc & aic/auc_brt_nest.csv")

###############################################################################
############################ GLM var & correlation ############################
###############################################################################

all_variables_30m <- c("aspect_30m", "brightness_30m", "cti_30m", "elev_30m", "evi_30m", "greenness_30m", "hli_30m", "roughness_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")

# correlation filter 0.75 / 0.8
glm_feed_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_hunt_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_nest_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")

# scale selection (AUC from uni-variate GLMs)
glm_feed_scale_all <- c("aspect_2000k", "cti_1000k", "elev_250k", "evi_1500k", "hli_1000k", "slope_250k", "tpi_2000k", "treecover_1500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_hunt_scale_all <- c("aspect_1000k", "cti_500k", "elev_30m", "evi_1500k", "hli_1500k", "slope_1000k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_nest_scale_all <- c("aspect_1500k", "cti_1000k", "elev_250k", "evi_1500k", "hli_100k", "slope_250k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")

# correlation filter 0.77 / 0.8 (-hli_1500k for hunt model)
glm_feed_scale <- c("aspect_2000k", "cti_1000k", "elev_250k", "evi_1500k", "hli_1000k", "slope_250k", "tpi_2000k", "treecover_1500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_hunt_scale <- c("aspect_1000k", "cti_500k", "elev_30m", "evi_1500k", "slope_1000k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
glm_nest_scale <- c("aspect_1500k", "cti_1000k", "elev_250k", "evi_1500k", "hli_100k", "slope_250k", "tpi_2000k", "treecover_1000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")

# Correlation plot V2
cor_list <- list()
for (i in 1:100) { # Loop through the 100th coord_feed elements
  extracted_data_feed <- list()
  for (variable in brt_feed_scale) {   # Extract data for variables
    raster_layer <- get(paste0(variable))
    extracted_data_feed[[variable]] <- extract(raster_layer, pa_feed_list[[i]][["coord_feed"]])
  }
  data_feed <- cbind(do.call(cbind, extracted_data_feed))   # Combine all extracted data into a dataframe
  env.cor <- cor(data_feed, method = "spearman")   # Compute the correlation matrix using Spearman method
  cor_list[[i]] <- env.cor   # Store the correlation matrix in the list
}
avg_cor_matrix_feed <- Reduce("+", cor_list) / length(cor_list) # Compute the average correlation matrix
corrplot(avg_cor_matrix_feed, method = "color", type = "upper", addCoef.col = 1, number.cex = 0.65, pch.cex = 1, diag = FALSE) # Plot the average correlation matrix

###############################################################################
################################### BRT var ################################### 
###############################################################################

brt_feed_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
brt_hunt_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
brt_nest_30m <- c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village")
brt_feed_scale <- c("aspect_1500k", "cti_1500k", "elev_30m", "evi_100k", "hli_2000k", "slope_250k", "tpi_2000k", "treecover_2000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
brt_hunt_scale <- c("aspect_2000k", "cti_1000k", "elev_30m", "evi_1500k", "hli_2000k", "slope_100k", "tpi_1000k", "treecover_500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village")
brt_nest_scale <- c("aspect_1500k", "cti_500k", "elev_100k", "evi_1000k", "hli_1500k", "slope_500k", "tpi_2000k", "treecover_2000k", "wetness_2000k", "distance_camp", "distance_river", "distance_road", "distance_village")

###############################################################################
################################# Dataframes ##################################
###############################################################################

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

df_glm_feed_30m_list <- create_data_frames_v2 (pa_feed_list, "feed", glm_feed_30m)
df_glm_hunt_30m_list <- create_data_frames_v2 (pa_hunt_list, "hunt", glm_hunt_30m)
df_glm_nest_30m_list <- create_data_frames_v2 (pa_nest_list, "nest", glm_nest_30m)

df_glm_feed_scale_list <- create_data_frames_v2 (pa_feed_list, "feed", glm_feed_scale)
df_glm_hunt_scale_list <- create_data_frames_v2 (pa_hunt_list, "hunt", glm_hunt_scale)
df_glm_nest_scale_list <- create_data_frames_v2 (pa_nest_list, "nest", glm_nest_scale)

df_brt_feed_30m_list <- create_data_frames_v2 (pa_feed_list, "feed", brt_feed_30m)
df_brt_hunt_30m_list <- create_data_frames_v2 (pa_hunt_list, "hunt", brt_hunt_30m)
df_brt_nest_30m_list <- create_data_frames_v2 (pa_nest_list, "nest", brt_nest_30m)

df_brt_feed_scale_list <- create_data_frames_v2 (pa_feed_list, "feed", brt_feed_scale)
df_brt_hunt_scale_list <- create_data_frames_v2 (pa_hunt_list, "hunt", brt_hunt_scale)
df_brt_nest_scale_list <- create_data_frames_v2 (pa_nest_list, "nest", brt_nest_scale)

# Define a named list of all the data frames to save
df_all_lists <- list(
  df_glm_feed_30m_list = df_glm_feed_30m_list,
  df_glm_hunt_30m_list = df_glm_hunt_30m_list,
  df_glm_nest_30m_list = df_glm_nest_30m_list,
  df_glm_feed_scale_list = df_glm_feed_scale_list,
  df_glm_hunt_scale_list = df_glm_hunt_scale_list,
  df_glm_nest_scale_list = df_glm_nest_scale_list,
  df_brt_feed_30m_list = df_brt_feed_30m_list,
  df_brt_hunt_30m_list = df_brt_hunt_30m_list,
  df_brt_nest_30m_list = df_brt_nest_30m_list,
  df_brt_feed_scale_list = df_brt_feed_scale_list,
  df_brt_hunt_scale_list = df_brt_hunt_scale_list,
  df_brt_nest_scale_list = df_brt_nest_scale_list
)

# Loop over the list and save each item
for (name in names(df_all_lists)) {
  save(list = name, file = paste0(name, ".RData"))
}

###############################################################################
############################## Dataframes 85-15 ###############################
###############################################################################

load("/Users/lucasdegreve/Data/RData/pa_list/pa_hunt_list85.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_nest_list85.RData")
load("/Users/lucasdegreve/Data/RData/pa_list/pa_feed_list85.RData")

df_glm_feed_30m_list85 <- create_data_frames_v2 (pa_feed_list85, "feed85", glm_feed_30m)
df_glm_hunt_30m_list85 <- create_data_frames_v2 (pa_hunt_list85, "hunt85", glm_hunt_30m)
df_glm_nest_30m_list85 <- create_data_frames_v2 (pa_nest_list85, "nest85", glm_nest_30m)

df_glm_feed_scale_list85 <- create_data_frames_v2 (pa_feed_list85, "feed85", glm_feed_scale)
df_glm_hunt_scale_list85 <- create_data_frames_v2 (pa_hunt_list85, "hunt85", glm_hunt_scale)
df_glm_nest_scale_list85 <- create_data_frames_v2 (pa_nest_list85, "nest85", glm_nest_scale)

df_85_lists <- list(
  df_glm_feed_30m_list85 = df_glm_feed_30m_list85,
  df_glm_hunt_30m_list85 = df_glm_hunt_30m_list85,
  df_glm_nest_30m_list85 = df_glm_nest_30m_list85,
  df_glm_feed_scale_list85 = df_glm_feed_scale_list85,
  df_glm_hunt_scale_list85 = df_glm_hunt_scale_list85,
  df_glm_nest_scale_list85 = df_glm_nest_scale_list85
)

# Loop over the list and save each item
for (name in names(df_85_lists)) {
  save(list = name, file = paste0(name, ".RData"))
}

###############################################################################
################################# GLM models ##################################
###############################################################################

library(MuMIn)
# Function to run dredge on a list of dataframes
run_dredge <- function(df_list) {
  dredge_results_list <- list()  # Initialize a list to store the dredge results for each dataframe
  for (i in 1:100) {  # Loop over the first 100 dataframes in df_list
    df <- df_list[[i]]  # Get the current dataframe
    formula <- as.formula(paste("presence ~", paste(names(df)[4:ncol(df)], collapse = " + ")))
    global_model <- glm(formula, data = df, family = binomial, na.action = na.fail)
    dredge_results <- dredge(global_model, trace = TRUE)
    dredge_results_list[[i]] <- dredge_results
  }
  return(dredge_results_list)
}
# List of dataset lists
df_glm_lists <- list(
  df_glm_feed_30m_list = df_glm_feed_30m_list,
  df_glm_hunt_30m_list = df_glm_hunt_30m_list,
  df_glm_nest_30m_list = df_glm_nest_30m_list,
  df_glm_feed_scale_list = df_glm_feed_scale_list,
  df_glm_hunt_scale_list = df_glm_hunt_scale_list,
  df_glm_nest_scale_list = df_glm_nest_scale_list
)
df_85_lists <- list(
  df_glm_feed_30m_list85 = df_glm_feed_30m_list85,
  df_glm_hunt_30m_list85 = df_glm_hunt_30m_list85,
  df_glm_nest_30m_list85 = df_glm_nest_30m_list85,
  df_glm_feed_scale_list85 = df_glm_feed_scale_list85,
  df_glm_hunt_scale_list85 = df_glm_hunt_scale_list85,
  df_glm_nest_scale_list85 = df_glm_nest_scale_list85
)

# df_glm_lists
for (name in names(df_85_lists)) {
  dredge_results_list <- run_dredge(df_85_lists[[name]])
  dredge_list_name <- paste0(name, "_dredge_results")
  assign(dredge_list_name, dredge_results_list)
  save(list = dredge_list_name, file = paste0(dredge_list_name, ".RData"))
  rm(list = dredge_list_name)
}

# Apply the function to each list and save the results
for (name in names(df_brt_lists)) {
  brt_results <- run_brt(df_brt_lists[[name]], elev_30m)
  variable_name <- paste0(name, "_brt_results")
  assign(variable_name, brt_results)
  save(list = variable_name, file = paste0(variable_name, ".RData"))
  rm(list = variable_name)
}

###############################################################################
################################# BRT models ##################################
###############################################################################

library(dismo)
library(blockCV)
n.folds <- 5
theRanges <- c(1300, 1300)

# Function to run BRT on a list of datasets
run_brt <- function(df_list, elev_30m) {
  brt_results <- list()   # Empty list to store results
  for (i in 1:100) {   # Loop over each dataset in the list
    dataset <- df_list[[i]]
    # Prepare spatial data
    spdf <- SpatialPointsDataFrame(coords = dataset[, c("x", "y")], data = dataset, proj4string = crs(elev_30m))
    # Generate cross-validation folds
    folds_with_similar_sizes <- FALSE
    while (!folds_with_similar_sizes) {
      myblocks <- cv_spatial(spdf, column = "presence", k = n.folds, size = theRanges, selection = "random", progress = FALSE)
      fold.vector3 <- myblocks$folds_ids
      fold.vector_presences <- fold.vector3[which(dataset[, "presence"] == 1)]
      counts <- hist(fold.vector_presences, plot = FALSE)$counts
      props <- counts[which(counts > 0)] / sum(counts)
      if (min(props) > 0.05) folds_with_similar_sizes <- TRUE
    }
    gbm.x = 4:length(dataset)
    # Fit BRT model
    brt_model <- gbm.step(dataset, gbm.x, gbm.y = 3, offset = NULL, fold.vector = fold.vector3,
                          tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.80,
                          site.weights = rep(1, dim(dataset)[1]), var.monotone = rep(0, length(gbm.x)), n.folds = 5,
                          prev.stratify = TRUE, family = "bernoulli", n.trees = 100, step.size = 10, max.trees = 10000,
                          tolerance.method = "auto", tolerance = 0.001, plot.main = TRUE, plot.folds = FALSE, verbose = TRUE,
                          silent = FALSE, keep.fold.models = FALSE, keep.fold.vector = FALSE, keep.fold.fit = FALSE)
    brt_results[[i]] <- brt_model
  }
  return(brt_results)
}

# List of dataset lists
df_brt_lists <- list(
  df_brt_feed_30m_list = df_brt_feed_30m_list,
  df_brt_hunt_30m_list = df_brt_hunt_30m_list,
  df_brt_nest_30m_list = df_brt_nest_30m_list,
  df_brt_feed_scale_list = df_brt_feed_scale_list,
  df_brt_hunt_scale_list = df_brt_hunt_scale_list,
  df_brt_nest_scale_list = df_brt_nest_scale_list
)

# Apply the function to each list and save the results
for (name in names(df_brt_lists)) {
  brt_results <- run_brt(df_brt_lists[[name]], elev_30m)
  variable_name <- paste0(name, "_brt_results")
  assign(variable_name, brt_results)
  save(list = variable_name, file = paste0(variable_name, ".RData"))
  rm(list = variable_name)
}

###############################################################################
############################### BRT prediction ################################
###############################################################################

library(gbm)
library(raster)

# Define the variable names for each BRT model type
brt_variable_names <- list(
  brt_feed_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_hunt_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_nest_30m = c("aspect_30m", "cti_30m", "elev_30m", "evi_30m", "hli_30m", "slope_30m", "tpi_30m", "treecover_30m", "wetness_30m", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_feed_scale = c("aspect_1500k", "cti_1500k", "elev_30m", "evi_100k", "hli_2000k", "slope_250k", "tpi_2000k", "treecover_2000k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_hunt_scale = c("aspect_2000k", "cti_1000k", "elev_30m", "evi_1500k", "hli_2000k", "slope_100k", "tpi_1000k", "treecover_500k", "wetness_1500k", "distance_camp", "distance_river", "distance_road", "distance_village"),
  brt_nest_scale = c("aspect_1500k", "cti_500k", "elev_100k", "evi_1000k", "hli_1500k", "slope_500k", "tpi_2000k", "treecover_2000k", "wetness_2000k", "distance_camp", "distance_river", "distance_road", "distance_village")
)
# Define the lists of rasters for each BRT model type
brt_rasters <- list(
  brt_feed_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  brt_hunt_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  brt_nest_30m = list(aspect_30m, cti_30m, elev_30m, evi_30m, hli_30m, slope_30m, tpi_30m, treecover_30m, wetness_30m, distance_camp, distance_river, distance_road, distance_village),
  brt_feed_scale = list(aspect_1500k, cti_1500k, elev_30m, evi_100k, hli_2000k, slope_250k, tpi_2000k, treecover_2000k, wetness_1500k, distance_camp, distance_river, distance_road, distance_village),
  brt_hunt_scale = list(aspect_2000k, cti_1000k, elev_30m, evi_1500k, hli_2000k, slope_100k, tpi_1000k, treecover_500k, wetness_1500k, distance_camp, distance_river, distance_road, distance_village),
  brt_nest_scale = list(aspect_1500k, cti_500k, elev_100k, evi_1000k, hli_1500k, slope_500k, tpi_2000k, treecover_2000k, wetness_2000k, distance_camp, distance_river, distance_road, distance_village)
)
# List of RData files to load
rdata_files <- list(
  df_brt_feed_30m_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_30m_list_brt_results.RData",
  df_brt_feed_scale_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_feed_scale_list_brt_results.RData",
  df_brt_hunt_30m_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_30m_list_brt_results.RData",
  df_brt_hunt_scale_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_hunt_scale_list_brt_results.RData",
  df_brt_nest_30m_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_30m_list_brt_results.RData",
  df_brt_nest_scale_list_brt_results = "/Users/lucasdegreve/Data/RData/brt_results/df_brt_nest_scale_list_brt_results.RData"
)

# Convert raster to points
centroids <- rasterToPoints(elev_30m, spatial = TRUE)
centroids_df <- data.frame(centroids@coords)

# Function to process each BRT model type
process_brt <- function(model_name, variable_names, rasters, rdata_file) {
  load(rdata_file)
  # Extract data from each raster layer and compile into a data frame
  extracted_data <- lapply(rasters, extract, y = centroids_df)
  data_raster <- cbind(centroids_df, do.call(cbind, extracted_data))
  colnames(data_raster) <- c("x", "y", variable_names)
  # Save the data raster
  save(data_raster, file = paste0("data_raster_", model_name, ".RData"))
  prediction_df <- data.frame(matrix(ncol = 100, nrow = nrow(data_raster)))
  for (i in 1:100) {
    best_model = get(ls()[grep(paste0("df_", model_name, "_list_brt_results"), ls())])[[i]]
    prediction = predict.gbm(best_model, data_raster, best_model$gbm.call$best.trees, type = "response", single.tree = FALSE)
    prediction_df[, i] <- prediction
  }
  save(prediction_df, file = paste0("prediction_df_", model_name, ".RData"))
}

# Apply the function to each model type
for (name in names(brt_variable_names)) {
  process_brt(name, brt_variable_names[[name]], brt_rasters[[name]], rdata_files[[paste0("df_", name, "_list_brt_results")]])
}

data_raster_brt_feed_30m <- data_raster
save(data_raster_brt_feed_30m, file = "data_raster_brt_feed_30m.RData")
rm(data_raster)

prediction_df_brt_feed_30m <- prediction_df
save(prediction_df_brt_feed_30m, file = "prediction_df_brt_feed_30m.RData")
rm(prediction_df)