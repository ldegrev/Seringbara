library(dismo)
library(blockCV)
n.folds <- 5
theRanges <- c(1300, 1300)

# Function to run BRT on a list of datasets
run_brt <- function(df_list, elev_30m) {
  brt_results <- vector("list", 100)  # Initialize list with 100 slots
  i <- 1
  while (i <= 100) {  # Continue until all 100 slots are filled with non-NULL values
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
    gbm.x <- 4:length(dataset)
    # Fit BRT model
    brt_model <- try(gbm.step(dataset, gbm.x, gbm.y = 3, offset = NULL, fold.vector = fold.vector3,
                              tree.complexity = 5, learning.rate = 0.001, bag.fraction = 0.80,
                              site.weights = rep(1, dim(dataset)[1]), var.monotone = rep(0, length(gbm.x)), n.folds = 5,
                              prev.stratify = TRUE, family = "bernoulli", n.trees = 100, step.size = 10, max.trees = 10000,
                              tolerance.method = "auto", tolerance = 0.001, plot.main = TRUE, plot.folds = FALSE, verbose = TRUE,
                              silent = FALSE, keep.fold.models = FALSE, keep.fold.vector = FALSE, keep.fold.fit = FALSE), silent = TRUE)
    
    if (inherits(brt_model, "try-error") || is.null(brt_model)) {
      cat("Model", i, "failed. Retrying...\n")
      next  # Skip the current iteration and try again
    }
    
    brt_results[[i]] <- brt_model
    i <- i + 1  # Move to the next index only if a valid model is created
  }
  
  return(brt_results)
}

# List of dataset lists
df_brt_lists <- list(
  df_brt_hunt_30m_list = df_brt_hunt_30m_list
)

# Apply the function to each list and save the results
for (name in names(df_brt_lists)) {
  brt_results <- run_brt(df_brt_lists[[name]], elev_30m)
  variable_name <- paste0(name, "_brt_results")
  assign(variable_name, brt_results)
  save(list = variable_name, file = paste0(variable_name, ".RData"))
  rm(list = variable_name)
}

df_brt_hunt_30m_list_brt_resultsV2 <- df_brt_hunt_30m_list_brt_results