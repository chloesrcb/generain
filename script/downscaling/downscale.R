muse <- FALSE

if (muse) {
  # Get the muse folder
  folder_muse <- "/home/serrec/work_rainstsimu/downscaling"
  path_to_python <- "/home/serrec/.pyenv/versions/3.9.18/bin/python3.9"
  setwd(folder_muse)
  # Load libraries and set theme
  source("utils.R")
  source("pinnEV.R")
  source("config.R")
} else {
  # Load libraries and set theme
  source("./script/load_libraries.R")
  source("./script/downscaling/pinnEV.R")
  path_to_python <- "/home/cserreco/.pyenv/versions/3.9.18/bin/python3.9"
}

library(sf)
library(sp)
library(geosphere)
library(reticulate)
library(lubridate)

py_version <- "3.9.18"
path <- paste0(reticulate::virtualenv_root(), "/pinnEV_env/bin/python")
Sys.setenv(RETICULATE_PYTHON =  path_to_python)
# Sys.setenv(RETICULATE_PYTHON = path) # Set Python interpreter
Sys.setenv(RETICULATE_LOG_LEVEL = "DEBUG")
tf_version = "2.13.1"
reticulate::use_virtualenv("pinnEV_env", required = T)

library(keras)
library(tensorflow)

#  DOWNSCALING WITH pinnEV

# Load table
output_file <- paste0(data_folder, "downscaling/downscaling_table.csv")
df <- read.csv(output_file, sep=";")
df$time <- as.POSIXct(df$time, tz="GMT")

# Add cyclical time features
df <- add_time_features(df, time_col = "time")

# Choose predictors
pred_cols <- c(paste0("X",1:27),
               "lon_Y","lat_Y","lon_X","lat_X",
               "hour_sin","hour_cos","minute_sin","minute_cos",
               "day_sin","day_cos","month_sin","month_cos")

# Build time and site reference
times <- sort(unique(df$time))
sites <- sort(unique(df$station))

# Y matrix (T x S)
Y_wide <- df %>%
  select(time, station, Y_obs) %>%
  pivot_wider(names_from = station, values_from = Y_obs) %>%
  arrange(time)

Y_mat <- as.matrix(Y_wide[, -1, drop=FALSE])
colnames(Y_mat) <- colnames(Y_wide)[-1]
stopifnot(nrow(Y_mat) == length(times))

# X array (T x S x d)
d <- length(pred_cols)
X_arr <- array(NA_real_, dim = c(length(times), length(sites), d),
               dimnames = list(NULL, sites, pred_cols))

for (j in seq_along(pred_cols)) {
  pj <- pred_cols[j]
  tmp <- df %>%
    select(time, station, !!sym(pj)) %>%
    pivot_wider(names_from = station, values_from = !!sym(pj)) %>%
    arrange(time)
  matj <- as.matrix(tmp[, -1, drop=FALSE])
  # Ensure same site order
  matj <- matj[, sites, drop=FALSE]
  X_arr[,,j] <- matj
}

# Split train/valid by year (last year = validation)
years <- year(times)
valid_year <- max(years, na.rm = TRUE)
train_t <- which(years < valid_year)
valid_t <- which(years == valid_year)

message("Train years: ", paste(sort(unique(years[train_t])), collapse=", "))
message("Valid year: ", valid_year)

# Standardize predictors using train only
sc <- standardize_3d_by_train(X_arr, train_t)
X_scaled <- sc$X_scaled

# ========== 2A) Fit EGPD (intensity part) ==========
# pinnEV can be trained on full Y with zeros (loss masks zeros),
# but we still must mask validation times with -1e10 for Y.train and vice versa.
Y.train <- Y.valid <- Y_mat
Y.train[valid_t, ] <- -1e10
Y.valid[train_t, ] <- -1e10

# Initial params: compute from your previous EGPD per-site fitting OR set reasonable defaults
# If you already computed mean_sigma/mean_kappa/mean_xi earlier, reuse them.
# Here: quick safe defaults (replace!)
if (!exists("mean_sigma")) mean_sigma <- 0.5
if (!exists("mean_kappa")) mean_kappa <- 0.4
if (!exists("mean_xi"))    mean_xi    <- 0.4

X.s <- list(X.nn.s = X_scaled)
X.k <- list(X.nn.k = X_scaled)

set.seed(1)
tf$random$set_seed(1L)

NN.fit <- eGPD.NN.train(
  Y.train = Y.train,
  Y.valid = Y.valid,
  X.s = X.s,
  X.k = X.k,
  type = "MLP",
  n.ep = 200,
  batch.size = 100,
  init.scale = mean_sigma,
  init.kappa = mean_kappa,
  init.xi = mean_xi,
  widths = c(6,3),
  seed = 1
)

pred <- eGPD.NN.predict(X.s = X.s, X.k = X.k, model = NN.fit$model)

# pred$pred.sigma etc should be (T x S)
stopifnot(all(dim(pred$pred.sigma) == dim(Y_mat)))
stopifnot(all(dim(pred$pred.kappa) == dim(Y_mat)))

# ========== 2B) Fit Occurrence model (baseline logistic) ==========
# Response: 1 if rain > 0 else 0
Y_occ <- ifelse(Y_mat > 0, 1, 0)

# Flatten X_scaled for a simple glm baseline:
# We create one row per (t,s). Use only train times for training.
X_flat <- matrix(X_scaled, nrow = dim(X_scaled)[1]*dim(X_scaled)[2], ncol = d, byrow = FALSE)
colnames(X_flat) <- pred_cols

# Map indices
Tn <- dim(X_scaled)[1]; Sn <- dim(X_scaled)[2]
time_idx <- rep(1:Tn, each = Sn)
site_idx <- rep(1:Sn, times = Tn)

occ_df <- data.frame(
  t = time_idx,
  s = site_idx,
  y = as.vector(Y_occ),
  X_flat
)

train_rows <- occ_df$t %in% train_t
valid_rows <- occ_df$t %in% valid_t

fit_logistic <- glm(y ~ . - t - s, data = occ_df[train_rows, ], family = binomial())
p_hat <- rep(NA_real_, nrow(occ_df))
p_hat[valid_rows] <- predict(fit_logistic, newdata = occ_df[valid_rows, ], type = "response")
p_hat[train_rows] <- predict(fit_logistic, newdata = occ_df[train_rows, ], type = "response")

# Reshape p_hat to (T x S)
P_hat <- matrix(p_hat, nrow = Tn, ncol = Sn, byrow = TRUE)

# ========== 2C) Generate downscaled rainfall ==========
# Hurdle generation:
# I ~ Bernoulli(P_hat), if I=0 => Y=0, else Y ~ eGPD(sigma_hat,kappa_hat,xi_hat)
set.seed(123)
I_sim <- matrix(rbinom(Tn*Sn, size=1, prob = as.vector(P_hat)),
                nrow = Tn, ncol = Sn, byrow = TRUE)

Y_sim <- matrix(0, nrow = Tn, ncol = Sn)
for (tt in 1:Tn) {
  for (ss in 1:Sn) {
    if (I_sim[tt, ss] == 1) {
      Y_sim[tt, ss] <- reGPD(
        n = 1,
        sigma = pred$pred.sigma[tt, ss],
        kappa = pred$pred.kappa[tt, ss],
        xi    = mean_xi  # keep xi fixed (or pred$pred.xi[tt,ss] if you trust it)
      )
    }
  }
}

# ========== 2D) Export results in long format ==========
df_out <- expand.grid(
  time = times,
  station = sites,
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
) %>%
  arrange(time, station) %>%
  mutate(
    p_occ = as.vector(P_hat),
    sigma = as.vector(pred$pred.sigma),
    kappa = as.vector(pred$pred.kappa),
    xi    = mean_xi,
    Y_sim = as.vector(Y_sim),
    Y_obs = as.vector(Y_mat)
  )

out_path <- paste0(data_folder, "downscaling/downscaled_results_long.csv")
write.csv(df_out, out_path, row.names = FALSE)
message("✅ downscaled results written to: ", out_path)

# Optional quick sanity checks
cat("\n--- sanity checks ---\n")
cat("Obs zero fraction:", mean(df_out$Y_obs == 0, na.rm=TRUE), "\n")
cat("Sim zero fraction:", mean(df_out$Y_sim == 0, na.rm=TRUE), "\n")
cat("Obs mean (pos only):", mean(df_out$Y_obs[df_out$Y_obs>0], na.rm=TRUE), "\n")
cat("Sim mean (pos only):", mean(df_out$Y_sim[df_out$Y_sim>0], na.rm=TRUE), "\n")

###############################################################################
# END
###############################################################################
