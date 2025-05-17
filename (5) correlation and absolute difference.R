library(raster)
library(terra)

coords <- centroids_MCP[, c("x", "y")]

# Initialize a list to store results
results <- list()

# Define the raster data pairs and their names
raster_pairs <- list(
  list("prediction_df_brt_feed_30m", "prediction_df_brt_hunt_30m"),
  list("prediction_df_brt_feed_scale", "prediction_df_brt_hunt_scale"),
  list("prediction_df_brt_nest_30m", "prediction_df_brt_hunt_30m"),
  list("prediction_df_brt_nest_scale", "prediction_df_brt_hunt_scale"),
  list("prediction_df_brt_nest_30m", "prediction_df_brt_feed_30m"),
  list("prediction_df_brt_nest_scale", "prediction_df_brt_feed_scale"),
  list("predictions_df_glm_feed_30m", "predictions_df_glm_hunt_30m"),
  list("predictions_df_glm_feed_scale", "predictions_df_glm_hunt_scale"),
  list("predictions_df_glm_nest_30m", "predictions_df_glm_hunt_30m"),
  list("predictions_df_glm_nest_scale", "predictions_df_glm_hunt_scale"),
  list("predictions_df_glm_nest_30m", "predictions_df_glm_feed_30m"),
  list("predictions_df_glm_nest_scale", "predictions_df_glm_feed_scale")
)

# Loop through each pair
for (pair in raster_pairs) {
  # Extract names
  name1 <- pair[[1]]
  name2 <- pair[[2]]
  pair_name <- paste(name1, "vs", name2, sep = "_")
  
  # Compute means for both raster datasets
  dataset1_means <- rowMeans(get(name1))
  dataset2_means <- rowMeans(get(name2))
  
  # Unlist if necessary
  if (is.list(dataset1_means)) dataset1_means <- unlist(dataset1_means)
  if (is.list(dataset2_means)) dataset2_means <- unlist(dataset2_means)
  
  # Create SpatRasters
  r1 <- terra::rast(cbind(coords, dataset1_means), crs = "EPSG:4326", type = "xyz")
  r2 <- terra::rast(cbind(coords, dataset2_means), crs = "EPSG:4326", type = "xyz")
  
  # Convert to RasterLayer
  r1_raster <- raster(r1)
  r2_raster <- raster(r2)
  
  # Stack and compute Pearson correlation
  b <- stack(r1_raster, r2_raster)
  a <- layerStats(b, "pearson", na.rm = TRUE)
  corr_value <- a$`pearson correlation coefficient`["dataset1_means", "dataset2_means"]
  
  # Compute absolute difference and its mean
  d <- abs(r1_raster - r2_raster)
  mean_abs_diff <- cellStats(d, "mean")
  
  # Store results
  results[[pair_name]] <- list(
    correlation_value = corr_value,
    mean_absolute_difference = mean_abs_diff
  )
}

# Print results
results

# Flatten the results into a data frame
results_df <- do.call(rbind, lapply(names(results), function(pair_name) {
  data.frame(
    Pair = pair_name,
    Correlation = results[[pair_name]]$correlation_value,
    Mean_Abs_Diff = results[[pair_name]]$mean_absolute_difference,
    stringsAsFactors = FALSE
  )
}))

# Export to a CSV file
write.csv(results_df, file = "results_summary.csv", row.names = FALSE)