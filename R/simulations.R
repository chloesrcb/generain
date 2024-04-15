# library(RandomFields)
# library(foreach)
# library(doParallel)

# generate_Z <- function(Z, i, modelBuhlCklu, Varm1, N, E, x, y, z) {
#   lx <- length(sx <- seq_along(x))  # spatial
#   ly <- length(sy <- seq_along(y))  # spatial
#   lz <- length(sz <- seq_along(z))  # temporal
#   # Z_i <- array(, dim = c(lx, ly, lz, 1))
#   Varm1_local <- Varm1
  
#   ## n=1
#   V <- 1/E[i,1] # processus de poisson
#   W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z, n=1) # processus gaussien
#   Y <- exp(W - W[1] - Varm1_local[,,,1])
#   Z[,,,i] <- V * Y

#   ## n dans {2,..,N}
#   for(n in 2:N) {
#     Exp <- E[i,n]
#     V <- 1/Exp
#     while(V > Z[N*(i-1)+n]) {
#       W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z)
#       Y <- exp(W - W[n] - Varm1_local[,,,n])
#       if(all(V*Y[seq_len(n-1)] < Z[(N*(i-1)+1):(N*(i-1)+(n-1))])) {
#           Z[,,,i] <- pmax(V*Y, Z[,,,i])
#       }
#       Exp <- Exp + rexp(1)
#       V <- 1/Exp
#     }
#   }
#   return(Z)
# }

# sim_BR_parallel <- function(beta1, beta2, alpha1, alpha2, x, y,
#                             z, n.BR, adv1 = 0, adv2 = 0) {
#   cl <- makeCluster(detectCores())
#   registerDoParallel(cl)

#   ## Setup
#   RFoptions(spConform = FALSE, install = "no", cores = detectCores())
#   lx <- length(sx <- seq_along(x))  # spatial
#   ly <- length(sy <- seq_along(y))  # spatial
#   lz <- length(sz <- seq_along(z))  # temporal
#   ## Variogram model 
#   modelBuhlCklu <- RMfbm(alpha = alpha1, var = beta1, proj = 1) +
#                    RMfbm(alpha = alpha1, var = beta1, proj = 2) +
#                    RMfbm(alpha = alpha2, var = beta2, proj = 3)

#   Nxy <- lx * ly
#   N <- Nxy * lz
#   grid <- matrix(0, nrow=N, ncol=3) # (N,3)-matrix

#   for (i in sx)
#     for (j in seq_len(ly*lz))
#       grid[i+(j-1)*ly, 1] <- i

#   for (i in sy)
#     for (j in sx)
#       for(k in sz)
#         grid[j + lx * (i-1) + (k-1) * Nxy, 2] <- i

#   for (i in sz)
#     for (j in seq_len(Nxy))
#       grid[j+Nxy*(i-1), 3] <- i

#   ## Shifted variogram
#   Varm1 <- foreach(n = seq_len(N), .combine = "c") %do% {
#     RandomFields::RFvariogram(modelBuhlCklu,
#                 x=(sx-grid[n,1]) - adv1*grid[n,3],
#                 y=(sy-grid[n,2]) - adv2*grid[n,3],
#                 z=sz-grid[n,3])
#   }
#   Varm1 <- array(Varm1, dim = c(lx, ly, lz, N))

#   #  Varm1 <- vapply(seq_len(N), function(n)
#   #               RandomFields::RFvariogram(modelBuhlCklu,
#   #                           x=(sx-grid[n,1]) - adv1*grid[n,3],
#   #                           y=(sy-grid[n,2]) - adv2*grid[n,3],
#   #                           z=sz-grid[n,3]),
#   #               array(NA_real_, dim=c(lx, ly, lz)))
#   ## => (lx, ly, lz, N)-array

#   ## Principal
#   Z <- array(, dim = c(lx, ly, lz, n.BR)) # 3d array
#   E <- matrix(rexp(n.BR * N), nrow=n.BR, ncol=N)

#   Z <- foreach(i = 1:n.BR, .combine = "c") %dopar% {
#     generate_Z(Z, i, modelBuhlCklu, Varm1, N, E, x, y, z)
#   }
#   stopCluster(cl)
#   return(Z)
# }
