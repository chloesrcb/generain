








test_that("get_criterion calculates the correct values", {
  # Create a sample dataframe
  df_result <- data.frame(
    beta1 = c(1, 2, 3),
    beta2 = c(4, 5, 6),
    alpha1 = c(7, 8, 9),
    alpha2 = c(10, 11, 12)
  )

  # Define true parameters
  true_param <- c(2, 5, 8, 11)

  # Call the function
  result <- get_criterion(df_result, true_param)

  # Expected results
  expected_mean <- c(mean(df_result$beta1), mean(df_result$beta2),
                     mean(df_result$alpha1), mean(df_result$alpha2))
  expected_rmse <- c(sqrt(mean((true_param[1] - df_result$beta1)^2)),
                    sqrt(mean((true_param[2] - df_result$beta2)^2)),
                     sqrt(mean((true_param[3] - df_result$alpha1)^2)),
                     sqrt(mean((true_param[4] - df_result$alpha2)^2)))
  expected_mae <- c(mean(abs(true_param[1] - df_result$beta1)),
                    mean(abs(true_param[2] - df_result$beta2)),
                    mean(abs(true_param[3] - df_result$alpha1)),
                    mean(abs(true_param[4] - df_result$alpha2)))

  # Compare results
  expect_equal(result$mean, expected_mean)
  expect_equal(result$rmse, expected_rmse)
  expect_equal(result$mae, expected_mae)
})






test_that("neg_ll without advection", {
  # # Create sample data
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:300

  # # Number of realizations
  # nres <- 10

  # # Folder
  # foldername <- paste0("./data/simulations_BR/sim_25s_300t/")

  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "br_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1.5, 1)
  # hmax <- sqrt(17)
  # tau_vect <- 0:10

  # df_lags <- get_lag_vectors(sites_coords, tau_max = 10)

  # q <- 0.8
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                                 type = "brownresnick")

  # nll <- neg_ll(true_param, simu_df, df_lags, quantile = q, excesses = excesses)

  # nll_hmax <- neg_ll(true_param, simu_df, df_lags, hmax = sqrt(17),
  #                   quantile = q, excesses = excesses)

  # expect_false(nll == nll_hmax)

  # # theorical likelihood
  # nmarg <- get_marginal_excess(simu_df, quantile = q)
  # pmarg <- nmarg / nrow(simu_df)
  # chi <- theoretical_chi(true_param, df_lags)

  # ll_df <- excesses
  # Tmax <- nrow(simu_df)
  # ll_df$Tobs <- Tmax - ll_df$tau # number of possible excesses
  # ll_df$chi <- chi$chi

  # ll_df$pchi <- 1 - pmarg * ll_df$chi

  # ll_df$non_excesses <- ll_df$Tobs - ll_df$kij # number of non-excesses
  # ll_df$ll <- ll_df$kij * log(ll_df$chi) +
  #             ll_df$non_excesses * log(ll_df$pchi)
  # ll_1 <- ll_df$kij[1] * log(ll_df$chi[1]) +
  #             ll_df$non_excesses[1] * log(ll_df$pchi[1])
  # expect_equal(ll_df$ll[1], ll_1)

  # nll_verif <- -sum(ll_df$ll, na.rm = TRUE)

  # expect_equal(nll, nll_verif)
})




test_that("neg_ll with advection", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:300

  # # Number of realizations
  # nres <- 10

  # # Folder
  # foldername <- paste0("./data/simulations_BR/sim_25s_300t_adv_0508/")

  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "br_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # adv <- c(0.5, 0.8)
  # true_param <- c(0.4, 0.2, 1.5, 1, adv)

  # hmax <- sqrt(17)
  # tau_vect <- 0:10

  # df_lags <- get_lag_vectors(sites_coords, tau_max = 10)

  # q <- 0.8
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                               type = "brownresnick")

  # nll <- neg_ll(params = true_param, data = simu_df, df_lags = df_lags,
  #                  quantile = q, excesses = excesses, hmax = hmax)

  # # theorical likelihood
  # nmarg <- get_marginal_excess(simu_df, quantile = q)
  # pmarg <- nmarg / nrow(simu_df)
  # excesses_hmax <- excesses[excesses$hnorm <= hmax, ]
  # chi <- theoretical_chi(true_param, excesses_hmax)

  # ll_df <- excesses_hmax
  # Tmax <- nrow(simu_df)
  # ll_df$chi <- chi$chi
  # ll_df$Tobs <- Tmax - ll_df$tau # number of possible excesses
  # ll_df$pchi <- 1 - pmarg * ll_df$chi

  # ll_df$non_excesses <- ll_df$Tobs - ll_df$kij # number of non-excesses
  # ll_df$ll <- ll_df$kij * log(ll_df$chi) +
  #             ll_df$non_excesses * log(ll_df$pchi)
  # ll_1 <- ll_df$kij[1] * log(ll_df$chi[1]) +
  #             ll_df$non_excesses[1] * log(ll_df$pchi[1])
  # expect_equal(ll_df$ll[1], ll_1)

  # nll_verif <- -sum(ll_df$ll, na.rm = TRUE)

  # expect_equal(nll, nll_verif)

  # # hmax == NA
  # hmax <- NA

  # df_lags <- get_lag_vectors(sites_coords, tau_max = 10,
  #                             hmax = hmax)

  # q <- 0.8
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                               type = "brownresnick")

  # nll <- neg_ll(params = true_param, data = simu_df, df_lags = df_lags,
  #                  quantile = q, excesses = excesses, hmax = hmax)

  # # theorical likelihood
  # nmarg <- get_marginal_excess(simu_df, quantile = q)
  # pmarg <- nmarg / nrow(simu_df)
  # chi <- theoretical_chi(true_param, df_lags)

  # ll_df <- excesses
  # Tmax <- nrow(simu_df)
  # ll_df$chi <- chi$chi
  # ll_df$Tobs <- Tmax - ll_df$tau # number of possible excesses
  # ll_df$pchi <- 1 - pmarg * ll_df$chi

  # ll_df$non_excesses <- ll_df$Tobs - ll_df$kij # number of non-excesses
  # ll_df$ll <- ll_df$kij * log(ll_df$chi) +
  #             ll_df$non_excesses * log(ll_df$pchi)
  # ll_1 <- ll_df$kij[1] * log(ll_df$chi[1]) +
  #             ll_df$non_excesses[1] * log(ll_df$pchi[1])
  # expect_equal(ll_df$ll[1], ll_1)

  # nll_verif <- -sum(ll_df$ll, na.rm = TRUE)

  # expect_equal(nll, nll_verif)
})



test_that("neg_ll with different advections", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:300

  # # Number of realizations
  # nres <- 10

  # # Folder
  # foldername <- paste0("./data/simulations_BR/sim_25s_300t/")

  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "br_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1.5, 1)

  # hmax <- NA
  # tau_vect <- 0:10

  # df_lags <- get_lag_vectors(sites_coords, tau_max = 10)

  # q <- 0.8
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                                   type = "brownresnick")

  # nll <- neg_ll(params = true_param, data = simu_df,
  #                  df_lags = df_lags,
  #                  quantile = q, excesses = excesses, hmax = hmax)

  # # with adv
  # adv <- c(0., 0.)

  # nll_adv0 <- neg_ll(params = c(true_param, adv), data = simu_df,
  #                  df_lags = df_lags,
  #                  quantile = q, excesses = excesses, hmax = hmax)

  # expect_equal(nll, nll_adv0)

  # adv <- c(0.8, 10)
  # nll_adv <- neg_ll(params = c(true_param, adv), data = simu_df,
  #                  df_lags = df_lags,
  #                  quantile = q, excesses = excesses, hmax = hmax)

  # expect_true(nll_adv > nll)

})


test_that("optim BR without advection", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:300

  # # Number of realizations
  # nres <- 10

  # # Folder
  # foldername <- paste0("./data/simulations_BR/sim_25s_300t/")

  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "br_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1.5, 1)
  # hmax <- NA
  # tau_vect <- 0:10

  # df_lags <- get_lag_vectors(sites_coords, tau_max = max(tau_vect))

  # q <- 0.95
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                                     type = "brownresnick")

  # result <- optim(par = true_param, fn = neg_ll,
  #               data = simu_df,
  #               quantile = q,
  #               df_lags = df_lags,
  #               excesses = excesses,
  #               hmax = hmax,
  #               method = "BFGS",
  #               control = list(parscale = c(1, 1, 1, 1),
  #                               maxit = 10000))
  # # check convergence
  # expect_true(result$convergence == 0)

  # # RMSE
  # df_result <- data.frame(
  #   beta1 = c(result$par[1]),
  #   beta2 = c(result$par[2]),
  #   alpha1 = c(result$par[3]),
  #   alpha2 = c(result$par[4])
  # )
  # df_valid <- get_criterion(df_result, true_param)

  # expect_true(df_valid$rmse[1] < 0.1)
  # expect_true(df_valid$rmse[2] < 0.1)
  # expect_true(df_valid$rmse[3] < 0.2)
  # expect_true(df_valid$rmse[4] < 0.2)

  # # adding null advection in optimization
  # result <- optim(par = c(true_param, 0., 0.), fn = neg_ll,
  #               data = simu_df,
  #               quantile = q,
  #               df_lags = df_lags,
  #               excesses = excesses,
  #               locations = sites_coords,
  #               hmax = hmax,
  #               method = "CG",
  #               control = list(parscale = c(1, 1, 1, 1, 1, 1),
  #                               maxit = 10000))
  # # check convergence
  # expect_true(result$convergence == 0)

  # # results
  # df_result <- data.frame(
  #   beta1 = c(result$par[1]),
  #   beta2 = c(result$par[2]),
  #   alpha1 = c(result$par[3]),
  #   alpha2 = c(result$par[4]),
  #   adv1 = c(result$par[5]),
  #   adv2 = c(result$par[6])
  # )
  # df_valid <- get_criterion(df_result, c(true_param, 0, 0))

  # expect_true(df_valid$rmse[1] < 0.1)
  # expect_true(df_valid$rmse[2] < 0.1)
  # expect_true(df_valid$rmse[3] < 0.25)
  # expect_true(df_valid$rmse[4] < 0.25)
  # expect_true(df_valid$rmse[5] < 0.25)
  # expect_true(df_valid$rmse[6] < 0.25)


  # # moving initial values
  # init_param <- c(true_param[1] + 0.1, true_param[2] + 0.2,
  #                 true_param[3], true_param[4])
  # result2 <- optim(par = init_param, fn = neg_ll,
  #               data = simu_df,
  #               quantile = q,
  #               df_lags = df_lags,
  #               excesses = excesses,
  #               locations = sites_coords,
  #               hmax = hmax,
  #               method = "BFGS",
  #               control = list(parscale = c(1, 1, 1, 1),
  #                               maxit = 10000))
  # # check convergence
  # expect_true(result$convergence == 0)

  # # results
  # df_result <- data.frame(
  #   beta1 = result2$par[1],
  #   beta2 = result2$par[2],
  #   alpha1 = result2$par[3],
  #   alpha2 = result2$par[4])
  
  # df_valid <- get_criterion(df_result, c(true_param))

  # expect_true(df_valid$rmse[1] < 0.1)
  # expect_true(df_valid$rmse[2] < 0.1)
  # expect_true(df_valid$rmse[3] < 0.25)
  # expect_true(df_valid$rmse[4] < 0.25)

  # # moving initial advection values
  # init_param <- c(true_param[1], true_param[2],
  #                 true_param[3], true_param[4], 0.1, 0.1)
  # df_lags <- get_lag_vectors(sites_coords, init_param, tau_vect = tau_vect)
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags)
  # result <- optim(par = init_param, fn = neg_ll,
  #               data = simu_df,
  #               quantile = q,
  #               df_lags = df_lags,
  #               excesses = excesses,
  #               locations = sites_coords,
  #               hmax = NA,
  #               method = "BFGS",
  #               control = list(parscale = c(1, 1, 1, 1, 1, 1),
  #                               maxit = 10000))
  # # check convergence
  # expect_true(result$convergence == 0)

  # # results
  # df_result <- data.frame(
  #   beta1 = c(result$par[1]),
  #   beta2 = c(result$par[2]),
  #   alpha1 = c(result$par[3]),
  #   alpha2 = c(result$par[4]),
  #   adv1 = c(result$par[5]),
  #   adv2 = c(result$par[6])
  # )
  # df_valid <- get_criterion(df_result, c(true_param, 0, 0))

  # expect_true(df_valid$rmse[1] < 0.1)
  # expect_true(df_valid$rmse[2] < 0.1)
  # expect_true(df_valid$rmse[3] < 0.25)
  # expect_true(df_valid$rmse[4] < 0.25)
  # expect_true(df_valid$rmse[5] < 0.25)
  # expect_true(df_valid$rmse[6] < 0.25)
})



test_that("optim BR with advection", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:300
  # adv <- c(0.5, 0.8)

  # # Number of realizations
  # nres <- 10
  # # Folder
  # foldername <- paste0("./data/simulations_BR/sim_25s_300t_adv_0105/")
  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "br_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1.5, 1, 0.1, 0.5)
  # hmax <- NA
  # tau_vect <- 0:8

  # df_lags <- get_lag_vectors(sites_coords, true_param, tau_vect = tau_vect)

  # q <- 0.9
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags)

  # result <- optim(par = true_param, fn = neg_ll,
  #               data = simu_df,
  #               quantile = q,
  #               df_lags = df_lags,
  #               excesses = excesses,
  #               locations = sites_coords,
  #               hmax = NA,
  #               method = "BFGS",
  #               control = list(parscale = c(1, 1, 1, 1, 1, 1),
  #                               maxit = 10000))
  # # check convergence
  # expect_true(result$convergence == 0)
})




test_that("optim rpareto", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:30
  # adv <- c(0.5, 0.3)
  # s0 <- c(1, 1)
  # t0 <- 1

  # # Number of realizations
  # nres <- 10
  # # Folder
  # foldername <- paste0("./data/simulations_rpar/sim_25s_30t_05_03/")
  # list_simu <- list()
  # for (i in 1:nres) {
  #   file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
  #                         length(temp), "t_", i, ".csv")
  #   list_simu[[i]] <- read.csv(file_name)
  # }

  # simu_df <- list_simu[[1]]

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1, 1.2, adv)
  # hmax <- NA
  # tau_vect <- 0:10

  # df_lags <- get_conditional_lag_vectors(sites_coords, true_param, s0, t0,
  #                                       tau_vect = tau_vect)
  # q <- 1
  # excesses <- empirical_excesses(simu_df, quantile = q, df_lags = df_lags,
  #                                threshold = TRUE)

  # # Combine all simulations
  # simu_all <- list_simu[[1]]
  # for (i in 2:nres){
  #   simu_all <- rbind(simu_all, list_simu[[i]])
  # }


  # df_lags <- get_conditional_lag_vectors(sites_coords, true_param, s0, t0,
  #                                       tau_vect = 0:10)

  # list_excesses <- list()
  # for (i in 1:nres) {
  #   excesses <- empirical_excesses(list_simu[[i]], quantile = q,
  #                                   df_lags = df_lags, threshold = TRUE)
  #   list_excesses[[i]] <- excesses
  # }

  # result <- optim(par = c(true_param), fn = neg_ll_composite,
  #                 list_simu = list_simu,
  #                 quantile = q,
  #                 df_lags = df_lags,
  #                 list_excesses = list_excesses,
  #                 locations = sites_coords,
  #                 hmax = 10,
  #                 s0 = s0,
  #                 t0 = t0,
  #                 threshold = TRUE,
  #                 method = "BFGS",
  #                 control = list(parscale = c(1, 1, 1, 1, 1, 1),
  #                                maxit = 10000))



})


test_that("optim on replicates with neg_ll_composite", {
  # ngrid <- 5
  # sites_coords <- generate_grid_coords(ngrid)
  # temp <- 1:30
  # adv <- c(0.5, 0.3)
  # s0 <- c(1, 1)
  # t0 <- 1

  # # Calculate the h vector
  # true_param <- c(0.4, 0.2, 1.5, 1, adv)
  # hmax <- NA
  # tau_vect <- 0:10

  # df_lags <- get_conditional_lag_vectors(sites_coords, true_param, s0, t0,
  #                                       tau_vect = tau_vect)

  # adv_int <- adv * 10
  # adv_str <- sprintf("%02d_%02d", adv_int[1], adv_int[2])

  # foldername <- paste0("../data/simulations_rpar/sim_", ngrid^2, "s_",
  #                             length(temp), "t_", adv_str, "/")

  # nres <- 10
  # list_simu <- list()
  # for (i in 1:nres) {
  #     file_name <- paste0(foldername, "rpar_", ngrid^2, "s_",
  #                             length(temp), "t_", i, ".csv")
  #     list_simu[[i]] <- read.csv(file_name)
  # }

  # # Combine all simulations
  # simu_all <- do.call(rbind, list_simu)

  # # simu_df_1 <- list_simu[[1]]
  # # simu_df_2 <- list_simu[[2]]
  # # simu_df_3 <- list_simu[[3]]

  # list_excesses <- list()
  # list_neg_ll <- list()
  # for (i in 1:nres) {
  #   list_excesses[[i]] <- empirical_excesses(list_simu[[i]], quantile = 1,
  #                                   df_lags = df_lags, threshold = TRUE)
  #   list_neg_ll[[i]] <- neg_ll(params = true_param, data = list_simu[[i]],
  #                               df_lags = df_lags, locations = sites_coords,
  #                               quantile = 1, excesses = list_excesses[[i]],
  #                               hmax = hmax, s0 = s0, t0 = t0, threshold = TRUE)
  #   print(list_neg_ll[[i]])
  # }

  # # neg_ll1 <- neg_ll(params = true_param, data = simu_df_1,
  # #                   df_lags = df_lags, locations = sites_coords,
  # #                   quantile = 1, excesses = list_excesses[[1]],
  # #                   hmax = hmax, s0 = s0, t0 = t0, threshold = TRUE)

  # # neg_ll2 <- neg_ll(params = true_param, data = simu_df_2,
  # #                   df_lags = df_lags, locations = sites_coords,
  # #                   quantile = 1, excesses = list_excesses[[2]],
  # #                   hmax = hmax, s0 = s0, t0 = t0, threshold = TRUE)

  # # neg_ll3 <- neg_ll(params = true_param, data = simu_df_3,
  # #                   df_lags = df_lags, locations = sites_coords,
  # #                   quantile = 1, excesses = list_excesses[[3]],
  # #                   hmax = hmax, s0 = s0, t0 = t0, threshold = TRUE)

  # neg_ll_all <- neg_ll_composite(params = true_param, list_simu = list_simu,
  #                               quantile = 1, df_lags = df_lags,
  #                               list_excesses = list_excesses,
  #                               locations = sites_coords, hmax = hmax,
  #                               s0 = s0, t0 = t0, threshold = TRUE)

  # # expect_equal(neg_ll_all, neg_ll1 + neg_ll2)

  # expect_equal(neg_ll_all, sum(unlist(list_neg_ll)))
})
