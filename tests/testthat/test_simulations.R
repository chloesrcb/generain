test_that("conditional_variogram", {
    x <- 1:10
    y <- 1:10
    t <- 1:30
    beta1 <- 0.8
    beta2 <- 0.2
    alpha1 <- 1.5
    alpha2 <- 1
    adv <- c(0, 0)

    lx <- length(sx <- seq_along(x))
    ly <- length(sy <- seq_along(y))
    lz <- length(sz <- seq_along(z))

    modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1,
                                                       proj = 1) +
                   RandomFields::RMfbm(alpha = alpha1, var = beta1,
                                       proj = 2) +
                   RandomFields::RMfbm(alpha = alpha2, var = beta2,
                                       proj = 3)

    ## Construct grid
    Nxy <- lx * ly
    N <- Nxy * lz
    grid <- matrix(0, nrow = N, ncol = 3) # (N,3)-matrix

    for (i in sx)
        for (j in seq_len(ly * lz))
            grid[i + (j - 1) * ly, 1] <- i

    for (i in sy)
        for (j in sx)
            for (k in sz)
                grid[j + lx * (i - 1) + (k - 1) * Nxy, 2] <- i

    for (i in sz)
        for (j in seq_len(Nxy))
            grid[j + Nxy * (i - 1), 3] <- i


    # Conditioning point
    s0 <- c(1, 1)
    s0_x <- s0[1]
    s0_y <- s0[2]
    t0 <- 1

    # Conditional semivariogram
    gamma_s0_t0 <- conditional_variogram(x, y, t, s0, t0, grid, modelBuhlCklu)

    expect_equal(dim(gamma_s0_t0), c(length(x), length(y), length(t)))

    # Check the first value
    expect_equal(gamma_s0_t0[s0_x, s0_y, t0], 0)

    # Check the second value
    s_x <- 1
    s_y <- 2
    time <- 1
    semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
                 beta1 * abs(s_y - s0_y)^alpha1 +
                 beta2 * abs(time - t0)^alpha2

    expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)

    # Check values
    s_x <- 2
    s_y <- 2
    time <- 2
    semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
                 beta1 * abs(s_y - s0_y)^alpha1 +
                 beta2 * abs(time - t0)^alpha2

    expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)

    # Check values
    s_x <- 3
    s_y <- 4
    time <- 8
    semivario_s_t <- beta1 * abs(s_x - s0_x)^alpha1 +
                 beta1 * abs(s_y - s0_y)^alpha1 +
                 beta2 * abs(time - t0)^alpha2

    expect_equal(gamma_s0_t0[s_x, s_y, time], semivario_s_t)


    # With advection
    adv <- c(0.5, 0.3)

    # check grid
    expect_equal(grid[1, 1], 1)
    expect_equal(grid[1, 2], 1)

    # Conditional semivariogram
    gamma_s0_t0_adv <- conditional_variogram(x, y, t, s0, t0,
                                         grid, modelBuhlCklu, adv)

    expect_equal(dim(gamma_s0_t0_adv), c(length(x), length(y), length(t)))

    # Check the first value
    semivario_s_t <- beta1 * abs(s0_x - s0_x - adv1 * 0)^alpha1 +
                 beta1 * abs(s0_y - s0_y - adv2 * 0)^alpha1 +
                 beta2 * abs(t0 - t0)^alpha2

    expect_equal(gamma_s0_t0[s0_x, s0_y, t0], semivario_s_t)

    # Check the other value
    s_x <- 1
    s_y <- 1
    time <- 2
    tau <- t0 - time
    semivario_s_t <- beta1 * abs(s_x - s0_x - adv[1] * tau)^alpha1 +
                 beta1 * abs(s_y - s0_y - adv[2] * tau)^alpha1 +
                 beta2 * abs(tau)^alpha2

    expect_equal(gamma_s0_t0_adv[s_x, s_y, time], semivario_s_t)

})

test_that("sim_rpareto", {
    # Test the sim_rpareto function
    x <- 1:10
    y <- 1:10
    t <- 1:30
    beta1 <- 0.4
    beta2 <- 0.2
    alpha1 <- 1.5
    alpha2 <- 1
    adv <- c(0.5, 0.3)

    W <- RandomFields::RFsimulate(modelBuhlCklu, x, y, t) # GP

    # conditional point
    s0 <- c(1, 1)
    t0 <- 1
    # check gaussian process
    expect_equal(dim(W), c(length(x), length(y), length(t)))
    # check conditional gaussian process W_s0,t0
    expect_equal(W[1, 1, t0], W[1]) # W_s0,t0 = W[s0[1], s0[2], t0]

    # Test the function
    Z <- sim_rpareto(beta1, beta2, alpha1, alpha2, x, y, t, adv, s0, t0, 10)

    # Check the length
    expect_equal(dim(Z), c(length(x), length(y), length(t), 10))

})


test_that("Check vario inside sim_BR", {
    # without advection
    x <- 1:10
    y <- 1:10
    z <- 1:30
    beta1 <- 0.4
    beta2 <- 0.2
    alpha1 <- 1.5
    alpha2 <- 1
    adv <- c(0, 0)

    lx <- length(sx <- seq_along(x))
    ly <- length(sy <- seq_along(y))
    lz <- length(sz <- seq_along(z))

    ## Model-Variogram BuhlCklu
    modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1,
                                                        proj = 1) +
                    RandomFields::RMfbm(alpha = alpha1, var = 2 * beta1,
                                        proj = 2) +
                    RandomFields::RMfbm(alpha = alpha2, var = 2 * beta2,
                                        proj = 3)

    ## Construct grid
    Nxy <- lx * ly
    N <- Nxy * lz
    grid <- matrix(0, nrow = N, ncol = 3) # (N,3)-matrix

    for (i in sx)
        for (j in seq_len(ly * lz))
            grid[i + (j - 1) * ly, 1] <- i

    for (i in sy)
        for (j in sx)
            for (k in sz)
                grid[j + lx * (i - 1) + (k - 1) * Nxy, 2] <- i

    for (i in sz)
        for (j in seq_len(Nxy))
            grid[j + Nxy * (i - 1), 3] <- i

    ## Construct shifted variogram
    gamma <- vapply(seq_len(N), function(n)
        RandomFields::RFvariogram(modelBuhlCklu,
            x = sx - grid[n, 1],
            y = sy - grid[n, 2],
            z = sz - grid[n, 3]),
            array(NA_real_, dim = c(lx, ly, lz))) ## => (lx, ly, lz, N)-array

    # check dimension
    expect_equal(dim(gamma), c(length(x), length(y), length(z), N))

    # Check value
    s1_x <- 1
    s1_y <- 1
    s2_x <- 1
    s2_y <- 1
    t1 <- 1
    t2 <- 1
    s2t2_index <- 1 + (s2_x - 1) + (s2_y - 1) * length(x) +
                  (t2 - 1) * length(x) * length(y)
    s1t1_index <- 1 + (s1_x - 1) + (s1_y - 1) * length(x) +
                  (t1 - 1) * length(x) * length(y)

    semivario_s_t <- beta1 * abs(s1_x - s2_x)^alpha1 +
                 beta1 * abs(s1_y - s2_y)^alpha1 +
                 beta2 * abs(t1 - t2)^alpha2

    expect_equal(gamma[s1_x,s1_y, t1, s2t2_index], 2 * semivario_s_t)
    expect_equal(gamma[s2_x,s2_y, t2, s1t1_index], 2 * semivario_s_t)

    s1_x <- 3
    s1_y <- 4
    s2_x <- 1
    s2_y <- 1
    t1 <- 10
    t2 <- 1
    s2t2_index <- 1 + (s2_x - 1) + (s2_y - 1) * length(x) +
                  (t2 - 1) * length(x) * length(y)
    s1t1_index <- 1 + (s1_x - 1) + (s1_y - 1) * length(x) +
                  (t1 - 1) * length(x) * length(y)

    semivario_s_t <- beta1 * abs(s1_x - s2_x)^alpha1 +
                 beta1 * abs(s1_y - s2_y)^alpha1 +
                 beta2 * abs(t1 - t2)^alpha2

    expect_equal(gamma[s1_x, s1_y, t1, s2t2_index], 2 * semivario_s_t)
    expect_equal(gamma[s2_x, s2_y, t2, s1t1_index], 2 * semivario_s_t)

    # With advection
    adv <- c(0.5, 0.3)

    # start.time <- Sys.time()
    # gamma <- array(0, dim = c(lx, ly, lz, N))

    # for (n in seq_len(N)) {
    #     for (i in seq_len(lx)) {
    #         for (j in seq_len(ly)) {
    #             for (k in seq_len(lz)) {
    #                 gamma[i, j, k, n] <- RandomFields::RFvariogram(
    #                     modelBuhlCklu,
    #                     x = x[i] - grid[n, 1] - adv[1] * (z[k] - grid[n, 3]),
    #                     y = y[j] - grid[n, 2] - adv[2] * (z[k] - grid[n, 3]),
    #                     z = z[k] - grid[n, 3]
    #                 )
    #                 }
    #             }
    #         }
    #     }
    # }
    # end.time <- Sys.time()
    # print(end.time - start.time)

    # # check dimension
    # expect_equal(dim(gamma), c(length(x), length(y), length(t), N))

    # # Check value
    # s1_x <- 1
    # s1_y <- 1
    # s2_x <- 1
    # s2_y <- 1
    # t1 <- 1
    # t2 <- 1
    # s2t2_index <- 1 + (s2_x - 1) + (s2_y - 1) * length(x) +
    #               (t2 - 1) * length(x) * length(y)
    # s1t1_index <- 1 + (s1_x - 1) + (s1_y - 1) * length(x) +
    #                 (t1 - 1) * length(x) * length(y)
    # semivario_s_t <- beta1 * abs(s1_x - s2_x - adv[1] * (t1 - t2))^alpha1 +
    #              beta1 * abs(s1_y - s2_y - adv[2] * (t1 - t2))^alpha1 +
    #              beta2 * abs(t1 - t2)^alpha2
    
    # expect_equal(gamma[s1_x, s1_y, t1, s2t2_index], semivario_s_t)


    # # Check value
    # s1_x <- 2
    # s1_y <- 3
    # s2_x <- 1
    # s2_y <- 2
    # t1 <- 10
    # t2 <- 1
    # s2t2_index <- 1 + (s2_x - 1) + (s2_y - 1) * length(x) +
    #               (t2 - 1) * length(x) * length(y)
    # s1t1_index <- 1 + (s1_x - 1) + (s1_y - 1) * length(x) +
    #                 (t1 - 1) * length(x) * length(y)
    # semivario_s_t <- beta1 * abs(s1_x - s2_x - adv[1] * (t1 - t2))^alpha1 +
    #              beta1 * abs(s1_y - s2_y - adv[2] * (t1 - t2))^alpha1 +
    #              beta2 * abs(t1 - t2)^alpha2

    # expect_equal(gamma[s1_x, s1_y, t1, s2t2_index], semivario_s_t)
})
