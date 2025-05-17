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
    
    # Determine y-axis limits
    y_min <- min(unlist(predictions_list))
    y_max <- max(unlist(predictions_list))
    
    plot(df[, covariates[i]], predictions_list[[1]], col="red", ann=FALSE, axes=FALSE, lwd=0.2, type="l", ylim = c(y_min, y_max))
    
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
