muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/downscaling"
  setwd(folder_muse)
  # Load libraries and set theme
  source("load_libraries.R")
  source("pinnEV.R")
  source("config.R")
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/downscaling/pinnEV.R")

}

library(sf)
library(sp)
library(geosphere)


# LOAD DOWNSCALING TABLE #######################################################

output_file <- paste0(data_folder, "downscaling/downscaling_table.csv")

df_Y_X_raw <- read.csv(output_file, sep = ";", header = TRUE)
head(df_Y_X_raw)
# colnames(df_Y_X_raw) <- c("time", "lon_Y", "lat_Y", "lon_X", "lat_X", "Y_obs",
#                           paste0("X", 1:27))

# Rename X_obs to X
# colnames(df_Y_X_raw)[which(colnames(df_Y_X_raw) == "X_obs")] <- "X"
nrow(df_Y_X_raw)
# head(df_Y_X_raw)

# first line
df_Y_X_raw[1,]
df_Y_X <- df_Y_X_raw

# Keep time as a POSIXct object
df_Y_X$time <- as.POSIXct(df_Y_X$time,
                         format = "%Y-%m-%d %H:%M:%S", tz = "GMT")

library(lubridate)
df_Y_X$hour <- hour(df_Y_X$time)
df_Y_X$minute <- minute(df_Y_X$time)
df_Y_X$day <- day(df_Y_X$time)
df_Y_X$month <- month(df_Y_X$time)
df_Y_X$year <- year(df_Y_X$time)
head(df_Y_X)


library(keras)
library(tensorflow)
library(reticulate)
py_version <- "3.9.18"
path <- paste0(reticulate::virtualenv_root(), "/pinnEV_env/bin/python")
Sys.setenv(RETICULATE_PYTHON = "/home/serrec/.virtualenvs/pinnEV_env/bin/python")
# Sys.setenv(RETICULATE_PYTHON = path) # Set Python interpreter
Sys.setenv(RETICULATE_LOG_LEVEL = "DEBUG")
tf_version = "2.13.1"
reticulate::use_virtualenv("pinnEV_env", required = T)

df_downscaling_all <- df_Y_X
# Remove rows with Y_obs = 0 for model EGPD of positive rainfall
df_downscaling_nn <- df_downscaling_all[df_downscaling_all$Y_obs > 0, ]
nrow(df_downscaling_nn)

n_observations <- nrow(df_downscaling_nn)
n_sites <- 17
n_predictors <- ncol(df_downscaling_nn) - 2

# remove the first column
df_downscaling_nn <- df_downscaling_nn[, -1]

# Standardize data
df_standardized <- df_downscaling_nn %>%
  mutate(across(where(is.numeric), scale))
head(df_standardized)

df_standardized <- df_downscaling_nn %>%
  # Convert time variables to cyclical features
  mutate(
    hour_sin = sin(2 * pi * hour / 24),
    hour_cos = cos(2 * pi * hour / 24),
    minute_sin = sin(2 * pi * minute / 60),
    minute_cos = cos(2 * pi * minute / 60),
    day_sin = sin(2 * pi * day / 31),  # Max 31 days
    day_cos = cos(2 * pi * day / 31),
    month_sin = sin(2 * pi * (month-1) / 12),
    month_cos = cos(2 * pi * (month-1) / 12)
  ) %>%
  select(-hour, -minute, -day, -month, -year) %>%  # Remove original columns
  # Standardize all numerical predictors (except cyclical time features)
  mutate(across(c(starts_with("X"), "Y_obs", "lon_X",
                        "lon_Y", "lat_X", "lat_Y"), scale))
head(df_standardized)
head(df_downscaling_nn)


X_all <- df_standardized[, -c(which(names(df_standardized) == "Y_obs"),
                           which(names(df_standardized) == "time"))]
# Convert into a 3D array
X_all <- array(unlist(X_all), dim = c(n_observations, n_sites, n_predictors))
Y_obs <- df_standardized$Y_obs

# Check for missing values
sum(is.na(X_all))
# sum(is.na(Y_obs))

# Identify rows with missing values in X_all
rows_with_na <- apply(X_all, 1, function(row) any(is.na(row)))

# Remove corresponding rows from X_all and Y_obs
X_all_clean <- X_all[!rows_with_na, , , drop = FALSE]  # Keep the same dimensions as the original array
Y_obs_clean <- Y_obs[!rows_with_na]
n_obs <- length(Y_obs_clean)
dim(X_all_clean)

# Flatten X_all_clean into a 2D array for glmnet
X_all_clean_flat <- array(X_all_clean, dim = c(n_obs, n_sites * n_predictors))
dim(X_all_clean_flat)

# Fit a LASSO model to select the most important predictors
library(glmnet)
fit_lasso <- cv.glmnet(as.matrix(X_all_clean_flat), Y_obs_clean, alpha = 1)


# Coefficients pour lambda.min (le meilleur lambda)
coef_lambda_min <- coef(fit_lasso, s = "lambda.min")
print(coef_lambda_min)
# Coefficients pour lambda.1se (le lambda à 1 écart-type)
coef(fit_lasso, s = "lambda.1se")

selected_predictors <- which(coef(fit_lasso, s = "lambda.min") != 0)
# Exclude intercept (index 1)
selected_predictors <- selected_predictors[selected_predictors > 1] - 1
X_all_clean_selected_flat <- X_all_clean_flat[, selected_predictors]


# Retransform the selected predictors into a 3D array
X_all_clean_selected <- array(X_all_clean_selected_flat, dim = c(n_obs, n_sites, length(selected_predictors)))
head(X_all_clean_selected)


# # Hyperparameter tuning for the neural network
# library(caret)

# grid <- expand.grid(
#   batch.size = c(10, 20),
#   widths = list(c(4, 2), c(8, 4, 2)),
#   n.ep = c(100, 200)
# )

# results <- apply(grid, 1, function(params) {
#   fit <- eGPD.NN.train(
#     Y.train, Y.valid, X.s, X.k, type = "MLP",
#     n.ep = params$n.ep, batch.size = params$batch.size,
#     widths = params$widths, seed = 1
#   )
#   return(fit$val_loss)  # Retourne la perte de validation
# })

# best_params <- grid[which.min(results), ]



X.nn.k <- X_all_clean # Nn predictors for kappa
X.nn.s <- X_all_clean # Nn predictors for sigma
length(Y_obs_clean)
sum(is.na(Y_obs_clean))

set.seed(123)
n_obs <- length(Y_obs_clean)
train_idx <- sample(seq_len(n_obs), size = 0.8 * n_obs)
Y.train <- Y_obs_clean
Y.valid <- Y_obs_clean
Y.train[-train_idx] <- -1e10
Y.valid[train_idx] <- -1e10
Y.train <- array(Y.train, dim = c(length(Y.train), 1))
Y.valid <- array(Y.valid, dim = c(length(Y.valid), 1))

# # NEW
# train_idx <- sample(seq_len(n_obs), size = 0.8 * n_obs)
# X.train <- X_all_clean_selected[train_idx, , ]
# X.valid <- X_all_clean_selected[-train_idx, , ]
# Y.train <- Y_obs_clean[train_idx]
# Y.valid <- Y_obs_clean[-train_idx]
# Y.train <- array(Y.train, dim = c(length(Y.train), 1))
# Y.valid <- array(Y.valid, dim = c(length(Y.valid), 1))

# X.nn.k <- X_all_clean_selected  # Nn predictors for kappa
# X.nn.s <- X_all_clean_selected # Nn predictors for sigma
# List of predictors
X.s <- list("X.nn.s" = X.nn.s)
X.k <- list("X.nn.k" = X.nn.k)


length(Y.train)

# initialize the neural network model EGPD fitting
# get censoring
censores <- seq(0.25, 0.35, 0.001)
rain_new <- rain_hsm[-1] # remove the first column date
# rain_cp <- rain_new[, c(5,6)]
df_score_ohsm <- choose_censore(rain_new, censores, n_samples = 100)

params_subrain <- get_egpd_estimates(rain_new,
                      left_censoring = 0.3)


mean_kappa <- mean(params_subrain$kappa)
mean_sigma <- mean(params_subrain$sigma)
mean_xi <- mean(params_subrain$xi)

# Train the neural network model
NN.fit <- eGPD.NN.train(Y.train, Y.valid, X.s, X.k, type = "MLP",
      n.ep = 200, batch.size = 10, init.scale = mean_sigma,
      init.kappa = mean_kappa, init.xi = mean_xi, widths = c(4, 2), seed = 1)
reticulate::py_last_error()
NN.fit
# we can change widths to reduce complexity (before c(6, 3))

# Predict with the trained model
out <- eGPD.NN.predict(X.s = X.s, X.k = X.k, NN.fit$model)
# out$pred.xi
# out$pred.kappa
# Results
print("sigma linear coefficients: "); print(round(out$pred.sigma, 5))
print("kappa linear coefficients: "); print(round(out$pred.kappa,2))

head(out$pred.sigma)

sigma_site1 <- out$pred.sigma[, 1]
kappa_site1 <- out$pred.kappa[, 1]

year_idx <- which(df_downscaling_nn$year == 2019)
sigma_site1_2019 <- sigma_site1[year_idx]
kappa_site1_2019 <- kappa_site1[year_idx]
print("Sigma for site 1 in 2019: ")
print(sigma_site1_2019)



################## OCCURENCE ###################################################
set.seed(123) # For reproducibility

df_downscaling_all <- df_Y_X
# Remove rows with Y_obs = 0 for model EGPD of positive rainfall
df_downscaling_nn <- df_downscaling_all[df_downscaling_all$Y_obs > 0, ]
nrow(df_downscaling_nn)

n_observations <- nrow(df_downscaling_nn)
n_sites <- 17
n_predictors <- ncol(df_downscaling_nn) - 2

# remove the first column
df_downscaling_nn <- df_downscaling_nn[, -1]


df_standardized <- df_downscaling_nn %>%
  # Convert time variables to cyclical features
  mutate(
    hour_sin = sin(2 * pi * hour / 24),
    hour_cos = cos(2 * pi * hour / 24),
    minute_sin = sin(2 * pi * minute / 60),
    minute_cos = cos(2 * pi * minute / 60),
    day_sin = sin(2 * pi * day / 31),  # Max 31 days
    day_cos = cos(2 * pi * day / 31),
    month_sin = sin(2 * pi * (month-1) / 12),
    month_cos = cos(2 * pi * (month-1) / 12)
  ) %>%
  select(-hour, -minute, -day, -month, -year) %>%  # Remove original columns
  # Standardize all numerical predictors (except cyclical time features)
  mutate(across(c(starts_with("X"), "lon_X",
                        "lon_Y", "lat_X", "lat_Y"), scale))
head(df_standardized)
head(df_downscaling_nn)


X_all <- df_standardized[, -c(which(names(df_standardized) == "Y_obs"),
                           which(names(df_standardized) == "time"))]
# Convert into a 3D array
X_all <- array(unlist(X_all), dim = c(n_observations, n_sites, n_predictors))
Y_obs <- df_standardized$Y_obs

# Check for missing values
sum(is.na(X_all))
# sum(is.na(Y_obs))

# Identify rows with missing values in X_all
rows_with_na <- apply(X_all, 1, function(row) any(is.na(row)))

# Remove corresponding rows from X_all and Y_obs
X_all_clean <- X_all[!rows_with_na, , , drop = FALSE]  # Keep the same dimensions as the original array
Y_obs_clean <- Y_obs[!rows_with_na]
n_obs <- length(Y_obs_clean)
dim(X_all_clean)

# Get occurrence of rainfall
Y_occurrence <- ifelse(Y_obs_clean > 0, 1, 0)

# Flatten dimensions 2 and 3 of X_all_clean into a data frame
X_train_df <- as.data.frame(matrix(X_all_clean[train_idx, , ], nrow = length(train_idx)))

# Rename columns
colnames(X_train_df) <- paste0("X", seq_len(ncol(X_train_df)))

# Create a data frame with the response variable and predictors
train_data <- data.frame(Y_occurrence = Y_occurrence[train_idx], X_train_df)

# Validation set
X_valid_df <- as.data.frame(matrix(X_all_clean[-train_idx, , ], nrow = length(Y_occurrence[-train_idx])))
colnames(X_valid_df) <- paste0("X", seq_len(ncol(X_valid_df)))

valid_data <- data.frame(Y_occurrence = Y_occurrence[-train_idx], X_valid_df)

# Logistic regression
fit_logistic <- glm(Y_occurrence ~ ., data = train_data, family = binomial)

summary(fit_logistic)

# Prediction on the validation set
valid_data$predicted_prob <- predict(fit_logistic, newdata = valid_data, type = "response")

# Convert the probabilities to classes (0 or 1) using a threshold of 0.5
valid_data$predicted_class <- ifelse(valid_data$predicted_prob > 0.5, 1, 0)

# # Evaluate the model
# library(caret)
# conf_matrix <- confusionMatrix(factor(valid_data$predicted_class), factor(valid_data$Y_occurrence))
# print(conf_matrix)

# ROC and AUC curves
# library(pROC)
# roc_curve <- roc(valid_data$Y_occurrence, valid_data$predicted_prob)
# plot(roc_curve, col = "blue", main = "ROC Curve")
# auc(roc_curve)
