# library(generain)
library(reshape2)
library(ggplot2)
setwd("./script")
# load libraries
source("load_libraries.R")



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


BR <- sim_BR(param[1], param[2], param[3], param[4], spa, spa, temp, n.BR)
save_simulations(BR, ngrid, n.BR,
        folder = paste0("../data/simulations_BR/oldsim_", ngrid^2, "s_",
                                length(temp), "t/"),
                  file = paste0("br_", ngrid^2, "s_",
                                length(temp), "t"), forcedind = 1)


file_path <- paste0("../data/simulations_BR/oldsim_25s_300t/br_",
                      ngrid^2, "s_", length(temp), "t_", 1, ".csv")
simu_df <- read.csv(file_path)

nsites <- ncol(simu_df)

sites_coords <- generate_grid_coords(sqrt(nsites))
dist_mat <- get_dist_mat(sites_coords,
                         latlon = FALSE) # distance matrix
df_dist <- reshape_distances(dist_mat) # reshape the distance matrix

true_param <- c(0.4, 0.2, 1.5, 1)
params <- true_param

sites_coords <- generate_grid_coords(sqrt(nsites))
h_vect <- get_lag_vectors(sites_coords, params,
                          hmax = sqrt(17), tau_vect = 0:10)
hmax <- sqrt(17)
q <- 0.9

chispa <- spatial_chi_alldist(df_dist, simu_df, quantile = q,
                                 hmax = hmax)
spa_estim <- get_estimate_variospa(chispa, weights = "exp", summary = FALSE)

q <- 0.9
tmax <- 10
chitemp <- temporal_chi(simu_df, tmax = tmax, quantile = q)
temp_params <- get_estimate_variotemp(chitemp, tmax, npoints = ncol(simu_df),
                                      weights = "exp", summary = FALSE)

df_result <- data.frame(beta1 =  spa_estim[1],
                        alpha1 = spa_estim[2],
                        beta2 = temp_estim[1],
                        alpha2 = temp_estim[2])
colnames(df_result) <- c("beta1", "alpha1", "beta2", "alpha2")

true_param <- c(beta1, beta2, alpha1, alpha2)
df_valid <- get_criterion(df_result, true_param)


# optim
q <- 0.9
excesses <- empirical_excesses(BR_df, q, h_vect)

true_param <- c(0.4, 0.2, 1.5, 1)
result <- optim(par = c(true_param), fn = neg_ll,
                        simu = BR_df,
                        quantile = q,
                        excesses = excesses,
                        h_vect = h_vect, tau = tau_vect,
                        locations = sites_coords,
                        method = "CG",
                        control = list(parscale = c(1, 1, 1, 1),
                                        maxit = 10000))
