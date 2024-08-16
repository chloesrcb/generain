library(generain)
library(reshape2)
library(ggplot2)
setwd("./script")
# load libraries
source("load_libraries.R")

# spatial and temporal structures
ngrid <- 5
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:300

# beta1, beta2, alpha1, alpha2
param <- c(0.8, 0.4, 1.5, 1) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
# beta3 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]
# alpha3 <- param[4]

adv <- c(0, 0)
true_param <- c(beta1, beta2, alpha1, alpha2)

# BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, 1,
#                adv)
# save_simulations(BR, ngrid, 1,
#                   folder = paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
#                                 length(temp), "t/"),
#                   file = paste0("br_", ngrid^2, "s_",
#                                 length(temp), "t"), forcedind = 1)


sim_BR <- function(beta1, beta2, alpha1, alpha2, x, y, z, n.BR) { 
  ## Setup 
  RandomFields::RFoptions(spConform=FALSE) 
  lx <- length(sx <- seq_along(x)) 
  ly <- length(sy <- seq_along(y)) 
  lz <- length(sz <- seq_along(z)) 
  ## Model-Variogram BuhlCklu 
  modelBuhlCklu <- RandomFields::RMfbm(alpha=alpha1, var=beta1, proj=1) + 
                  RandomFields::RMfbm(alpha=alpha1, var=beta1, proj=2) + 
                  RandomFields::RMfbm(alpha=alpha2, var=beta2, proj=3)
  
  ## Construct grid 
  Nxy <- lx * ly 
  N <- Nxy * lz 
  grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix 

  for (i in sx) 
    for (j in seq_len(ly*lz)) 
      grid[i+(j-1)*ly, 1] <- i 
  
  for (i in sy) 
    for (j in sx) 
      for(k in sz) 
        grid[j+lx*(i-1)+(k-1)*Nxy, 2] <- i 
  
  for (i in sz) 
    for (j in seq_len(Nxy)) 
      grid[j+Nxy*(i-1), 3] <- i

  ## Construct shifted variogram
  Varm1 <- vapply(seq_len(N), function(n) 
      RandomFields::RFvariogram(modelBuhlCklu, 
        x=sx-grid[n,1], 
        y=sy-grid[n,2], 
        z=sz-grid[n,3]), 
        array(NA_real_, dim=c(lx, ly, lz))) ## => (lx, ly, lz, N)-array 
  
  ## Main 
  Z <- array(, dim=c(lx, ly, lz, n.BR)) # 4d array 
  E <- matrix(rexp(n.BR * N), nrow=n.BR, ncol=N) 
  for (i in seq_len(n.BR)) { ## n=1 
    V <- 1/E[i,1] 
    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n=1) 
    Y <- exp(W - W[1] - Varm1[,,,1]) 
    Z[,,,i] <- V * Y 
    ## n in {2,..,N} 
    for(n in 2:N) { 
      Exp <- E[i,n] 
      V <- 1/Exp 
      while(V > Z[N*(i-1)+n]) { 
        W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z) 
        Y <- exp(W - W[n] - Varm1[,,,n]) 
        if(all(V*Y[seq_len(n-1)] < Z[(N*(i-1)+1):(N*(i-1)+(n-1))])) 
          Z[,,,i] <- pmax(V*Y, Z[,,,i]) 
          Exp <- Exp + rexp(1) 
          V <- 1/Exp 
      } 
    } 
  } 
  ## Return 
  Z 
}

library(foreach)
library(doParallel)

n.BR <- 1
num_iterations <- 5

# Create a parallel backend
cl <- makeCluster(detectCores())
registerDoParallel(cl)

# Run the simulation in parallel
results <- foreach(i = 1:num_iterations, .combine = "c") %dopar% {
  BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR)
  save_simulations(BR, ngrid, n.BR,
                  folder = paste0("../data/simulations_BR/oldsim_", ngrid^2, "s_",
                                length(temp), "t/"),
                  file = paste0("br_", ngrid^2, "s_",
                                length(temp), "t"), forcedind = i)

}

# Stop the parallel backend
stopCluster(cl)


# plot(BR[5,4,,1], main = "BR simulation")


# save simulations to CSV files
# save_simulations(BR, ngrid, n.BR,
#                  folder = paste0("../data/simulations_BR/sim_", ngrid^2, "s_",
#                                 length(temp), "t/"),
#                  file = paste0("br_", ngrid^2, "s_",
#                                 length(temp), "t"))


# Simulation
num_iterations <- 1
list_BR <- list()
for (i in 1:num_iterations) {
  file_path <- paste0("../data/simulations_BR/oldsim_25s_300t/br_",
                      ngrid^2, "s_", length(temp), "t_", i, ".csv")
  df <- read.csv(file_path)
  list_BR[[i]] <- df
}

# dependece buhl
# simu_df <- BR
simu_df <- list_BR[[1]] # first simulation
# simu_df <- simulation_data
par(mfrow=c(1,1))
nsites <- ncol(simu_df) # number of sites
plot(simu_df[, 1])
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix


params <- true_param
nsites <- ncol(simu_df) # number of sites
# get grid coordinates
sites_coords <- generate_grid_coords(sqrt(nsites))
# sites_coords$id <- 1:nsites
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 0)

# Evaluate the estimates with Buhl method
spa_estim <- evaluate_vario_estimates(list_BR, 0.8,
                                  spatial = TRUE, df_dist = df_dist,
                                  hmax = sqrt(17))

temp_estim <- evaluate_vario_estimates(list_BR, 0.85,
                                          spatial = FALSE, tmax = 10)

df_result <- data.frame(beta1 =  unlist(spa_estim$beta),
                        alpha1 = unlist(spa_estim$alpha),
                        beta2 = unlist(temp_estim$beta),
                        alpha2 = unlist(temp_estim$alpha))
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

true_param <- c(beta1, beta2, alpha1, alpha2)
df_valid <- get_criterion(df_result, true_param)

# save the results in latex table
write.table(round(df_valid, 3), file = "valid_br.tex",
            row.names = c("Beta1", "Alpha1", "Beta2", "Alpha2"),
            col.names = c("Mean", "RMSE", "MAE"),
            quote = FALSE)


################################################################################
# VERIFICATION
################################################################################


# qq-plots margins
par(mfrow = c(2, 2))


library(ismev)

BR_df <- simulation_data
BR_loc <- simu_df$S5
plot(BR_loc)
BR_loc_log <- log(BR_loc)
gumbel.fit <- gum.fit(BR_loc_log)

mu <- gumbel.fit$mle[1]
sigma <- gumbel.fit$mle[2]

library(evd)

theorical_qgum <- qgumbel(ppoints(BR_loc_log), mu, sigma)

qqplot(BR_loc_log, theorical_qgum, main = "Gumbel Q-Q plot",
       xlab = "Empirical quantiles",
       ylab = "Theoretical quantiles")


# fit frechet
library(extRemes)
frechet.fit <- fevd(BR_loc, type = "Gumbel", method = "MLE")

mu <- frechet.fit$results$par[1]
sigma <- frechet.fit$results$par[2]

theorical_qfrechet <- qgev(ppoints(BR_loc), mu, sigma, 0)

qqplot(BR_loc, theorical_qfrechet, main = "Frechet Q-Q plot",
  xlab = "Empirical quantiles",
  ylab = "Theoretical quantiles")
