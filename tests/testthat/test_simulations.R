# test_that("conditional_variogram", {
#     x <- 1:10
#     y <- 1:10
#     z <- 1:30
#     beta1 <- 0.8
#     beta2 <- 0.2
#     alpha1 <- 1.5
#     alpha2 <- 1
#     adv <- c(0, 0)

#     lx <- length(sx <- seq_along(x))
#     ly <- length(sy <- seq_along(y))
#     lz <- length(sz <- seq_along(z))

#     modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1,
#                                                        proj = 1) +
#                    RandomFields::RMfbm(alpha = alpha1, var = beta1,
#                                        proj = 2) +
#                    RandomFields::RMfbm(alpha = alpha2, var = beta2,
#                                        proj = 3)

#     ## Construct grid
#     Nxy <- lx * ly
#     N <- Nxy * lz
#     grid <- matrix(0, nrow = N, ncol = 3) # (N,3)-matrix

#     for (i in sx)
#         for (j in seq_len(ly * lz))
#             grid[i + (j - 1) * ly, 1] <- i

#     for (i in sy)
#         for (j in sx)
#             for (k in sz)
#                 grid[j + lx * (i - 1) + (k - 1) * Nxy, 2] <- i

#     for (i in sz)
#         for (j in seq_len(Nxy))
#             grid[j + Nxy * (i - 1), 3] <- i


#     # Conditioning point
#     s0 <- c(1, 1)
#     s0_x <- s0[1]
#     s0_y <- s0[2]
#     t0 <- 10
#     index_s0_t0 <- 1 + (s0_x - 1) + (s0_y - 1) * length(x) +
#                    (t0 - 1) * length(x) * length(y)

#     # Conditional semivariogram
#     gamma_s0_t0 <- conditional_variogram(x, y, z, s0, t0, grid, modelBuhlCklu)

#     expect_equal(dim(gamma_s0_t0), c(length(x), length(y), length(z)))

#     semivario_s_t <- beta1 * abs(s0_x - s0_x)^alpha1 +
#                  beta1 * abs(s0_y - s0_y)^alpha1 +
#                  beta2 * abs(t0 - t0)^alpha2

#     gamma <- vapply(seq_len(N), function(n)
#         RandomFields::RFvariogram(modelBuhlCklu,
#             x = sx - grid[n, 1],
#             y = sy - grid[n, 2],
#             z = sz - grid[n, 3]),
#             array(NA_real_, dim = c(lx, ly, lz))) ## => (lx, ly, lz, N)-array

#     gamma_s0_t0_bis <- gamma[, , , index_s0_t0]

#     # Check the first value
#     expect_equal(gamma_s0_t0[s0_x, s0_y, t0], semivario_s_t)
#     expect_equal(gamma_s0_t0[s0_x, s0_y, t0], gamma_s0_t0_bis[s0_x, s0_y, t0])

#     # Check the second value
#     s_x <- 1
#     s_y <- 2
#     time <- 1
#     semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
#                  beta1 * abs(s_y - s0_y)^alpha1 +
#                  beta2 * abs(time - t0)^alpha2

#     expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)

#     # Check values
#     s_x <- 2
#     s_y <- 2
#     time <- 2
#     semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
#                  beta1 * abs(s_y - s0_y)^alpha1 +
#                  beta2 * abs(time - t0)^alpha2

#     expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)

#     # Check values
#     s_x <- 3
#     s_y <- 4
#     time <- 8
#     semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
#                  beta1 * abs(s_y - s0_y)^alpha1 +
#                  beta2 * abs(time - t0)^alpha2

#     expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)


#     # With advection
#     adv <- c(0.5, 0.3)

#     # check grid
#     expect_equal(grid[1, 1], 1)
#     expect_equal(grid[1, 2], 1)

#     # Conditional semivariogram
#     gamma_s0_t0_adv <- conditional_variogram(x, y, z, s0, t0,
#                                          grid, modelBuhlCklu, adv)

#     expect_equal(dim(gamma_s0_t0_adv), c(length(x), length(y), length(z)))

#     # Check the first value
#     semivario_s_t <- beta1 * abs(s0_x - s0_x - adv[1] * 0)^alpha1 +
#                  beta1 * abs(s0_y - s0_y - adv[2] * 0)^alpha1 +
#                  beta2 * abs(t0 - t0)^alpha2

#     expect_equal(gamma_s0_t0[s0_x, s0_y, t0], semivario_s_t)

#     # Check the other value
#     s_x <- 1
#     s_y <- 1
#     time <- 2
#     tau <- t0 - time
#     semivario_s_t <- beta1 * abs(s_x - s0_x - adv[1] * tau)^alpha1 +
#                  beta1 * abs(s_y - s0_y - adv[2] * tau)^alpha1 +
#                  beta2 * abs(tau)^alpha2

#     expect_equal(gamma_s0_t0_adv[s_x, s_y, time], semivario_s_t)

# })

# test_that("sim_rpareto", {
#     # Test the sim_rpareto function
#     x <- 1:10
#     y <- 1:10
#     t <- 1:30
#     beta1 <- 0.4
#     beta2 <- 0.2
#     alpha1 <- 1.5
#     alpha2 <- 1
#     adv <- c(0.5, 0.3)

#     modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1,
#                                                        proj = 1) +
#                    RandomFields::RMfbm(alpha = alpha1, var = beta1,
#                                        proj = 2) +
#                    RandomFields::RMfbm(alpha = alpha2, var = beta2,
#                                        proj = 3)

#     # simulate the gaussian process
#     # W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, z=t)

#     # # conditional point
#     # s0 <- c(1, 1)
#     # t0 <- 1

#     # # check gaussian process
#     # expect_equal(dim(W), c(length(x), length(y), length(t)))
#     # # check conditional gaussian process W_s0,t0
#     # expect_equal(W[s0[1], s0[2], t0], W[1]) # W_s0,t0 = W[s0[1], s0[2], t0]

#     # Z <- sim_rpareto(beta1, beta2, alpha1, alpha2, x, y, t, adv, nres=10)
#     # # Check the length
#     # expect_equal(dim(Z), c(length(x), length(y), length(t), 10))

#     # # conditional point
#     # s0 <- c(2, 3)
#     # t0 <- 3
#     # index_s0_t0 <- 1 + (s0[1] - 1) + (s0[2] - 1) * length(x) +
#     #                (t0 - 1) * length(x) * length(y)

#     # # check conditional gaussian process W_s0,t0
#     # expect_equal(W[s0[1], s0[2], t0], W[index_s0_t0])

# })


# test_that("concat r-pareto simulations", {
#     # # Test the sim_rpareto function
#     # x <- 1:10
#     # y <- 1:10
#     # t <- 1:30
#     # beta1 <- 0.4
#     # beta2 <- 0.2
#     # alpha1 <- 1.5
#     # alpha2 <- 1
#     # adv <- c(0.5, 0.3)

#     # adv_int <- adv * 10
#     # adv_str <- sprintf("%02d_%02d", adv_int[1], adv_int[2])

#     # foldername <- paste0("../data/simulations_rpar/sim_", ngrid^2, "s_",
#     #                             length(temp), "t_", adv_str, "/")

#     # nres <- 2
#     # list_simu <- list()
#     # for (i in 1:nres) {
#     #     file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
#     #                             length(temp), "t_", i, ".csv")
#     #     list_simu[[i]] <- read.csv(file_name)
#     # }

#     # # Combine all simulations
#     # simu_all <- do.call(rbind, list_simu)

#     # simu_df_1 <- list_simu[[1]]
#     # nb_excesses_simu1 <- sum(simu_df_1$S1 > 1)
#     # simu_df_2 <- list_simu[[2]]
#     # nb_excesses_simu2 <- sum(simu_df_2$S1 > 1)
#     # nb_exccesses_all <- sum(simu_all$S1 > 1)

#     # expect_equal(nb_excesses_simu1 + nb_excesses_simu2, nb_exccesses_all)
# })
