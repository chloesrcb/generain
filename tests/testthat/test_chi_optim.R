# Description: This file contains the tests for the optimization functions.

# Example parameters and data for testing
eta1 <- 0.8
eta2 <- 0.3

df_lags <- data.frame(
  s1x = c(0, 1, 2, 0),
  s1y = c(0, 1, 2, 0),
  s2x = c(0, 2, 3, 0),
  s2y = c(0, 2, 3, 0),
  tau = c(0, 1, 2, 1),
  s1 = c(1, 2, 3, 1),
  s2 = c(1, 3, 4, 1),
  hnorm = c(0, 1, 2, 0)
)


### Tests for the theoretical_chi function -------------------------------------
test_that("Test that it works well for different configurations", {
  wind_vect <- c(10, 5)  # Example wind vector (vx = 10 m/s, vy = 5 m/s)
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  params <- c(0.1, 0.4, 0.5, 1.5, adv)  # Example variogram parameters
  # Test if the function returns a data frame
  result <- theoretical_chi(params, df_lags, latlon = FALSE, distance = "lalpha")
  expect_true(is.data.frame(result)) # Check if the result is a data frame
  # Check if the columns are present
  expect_true("chi" %in% colnames(result))
  expect_true("hlag" %in% colnames(result))
  expect_true("vario" %in% colnames(result))
  expect_true("hnormV" %in% colnames(result))

  # Test that the 'hnormV' values are positive or zero
  expect_true(all(result$hnormV >= 0))
  # Test that chi values are between 0 and 1 (as they are probabilities)
  expect_true(all(result$chi >= 0 & result$chi <= 1))
  # Test that the 'vario' values are positive or zero
  expect_true(all(result$vario >= 0))

  # Test that there are no NA or NaN values in critical columns
  expect_false(any(is.na(result$hlag)))
  expect_false(any(is.na(result$chi)))
  expect_false(any(is.na(result$vario)))
  expect_false(any(is.na(result$hlag)))
  expect_true(all(!is.na(result$x_polar)))
  expect_true(all(!is.na(result$y_polar)))

  # For same site pairs, the hnormV should be zero for tau = 0
  expect_true(all(result$hnormV[result$s1 == result$s2 & result$tau == 0] == 0))

  # For same site with different tau, the hnormV should be different
  expect_true(all(result$hnormV[result$s1 == result$s2 & result$tau != 0] != 0))

  # Optionally, test the behavior with latlon = TRUE (if needed)
  result_latlon <- theoretical_chi(params, df_lags,
                                    latlon = TRUE, distance = "lalpha")
  expect_true("chi" %in% colnames(result_latlon))
  expect_true("hlag" %in% colnames(result_latlon))
  expect_true("vario" %in% colnames(result_latlon))
  expect_true("hnormV" %in% colnames(result_latlon))
  expect_true("x_polar" %in% colnames(result_latlon))
  expect_true("y_polar" %in% colnames(result_latlon))

  # Test that the 'hnormV' values are positive or zero
  expect_true(all(result_latlon$hnormV >= 0))
  # Test that chi values are between 0 and 1 (as they are probabilities)
  expect_true(all(result_latlon$chi >= 0 & result$chi <= 1))
  # Test that the 'vario' values are positive or zero
  expect_true(all(result_latlon$vario >= 0))

  # Test that there are no NA or NaN values in critical columns
  expect_false(any(is.na(result_latlon$hlag)))
  expect_false(any(is.na(result_latlon$chi)))
  expect_false(any(is.na(result_latlon$vario)))
  expect_false(any(is.na(result_latlon$hlag)))
  expect_true(all(!is.na(result_latlon$x_polar)))
  expect_true(all(!is.na(result_latlon$y_polar)))

  result_nodir <- theoretical_chi(params, df_lags,
                                  latlon = FALSE, distance = "euclidean")
  expect_false("x_polar" %in% colnames(result_nodir))
  expect_false("y_polar" %in% colnames(result_nodir))

  # If no wind
  wind_vect <- c(0, 0)  # Example wind vector (vx = 10 m/s, vy = 5 m/s)
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  params <- c(0.1, 0.4, 0.5, 1.5, adv)  # Example variogram parameters
  result_nowind <- theoretical_chi(params, df_lags, latlon = FALSE,
                                   distance = "lalpha")

  # For same site pairs, the hnormV should be zero for tau = 0
  expect_true(all(result_nowind$hnormV[result_nowind$s1 == result_nowind$s2 & 
                                                  result_nowind$tau == 0] == 0))
  # For same site with different tau, the hnormV should be different than
  # the one with wind
  expect_true(result_nowind$hnormV[2] != result$hnormV[2])

  # If no wind vector is provided and no advection
  adv <- c(0, 0)
  params[5:6] <- adv
  result_nowind <- theoretical_chi(params, df_lags, latlon = FALSE,
                                   distance = "lalpha")

  # For same site pairs, the hnormV should be zero for tau = 0
  expect_true(all(result_nowind$hnormV[result_nowind$s1 == result_nowind$s2 &
                                                  result_nowind$tau == 0] == 0))
  # For same site with different tau, the hnormV should be zero bc no adv
  expect_true(all(result_nowind$hnormV[result_nowind$s1 == result_nowind$s2 &
                                                  result_nowind$tau == 1] == 0))

})


test_that("theoretical_chi returns correct structure", {
  params <- c(0.5, 0.3, 1.2, 0.8, 0.2, 0.1)
  df_lags <- data.frame(
    s1 = c(1, 2),
    s2 = c(2, 3),
    tau = c(1, 2),
    s1x = c(0, 1),
    s1y = c(0, 1),
    s2x = c(1, 2),
    s2y = c(1, 2),
    hnorm = c(1.41, 1.41)
  )

  result <- theoretical_chi(params, df_lags, latlon = FALSE, distance = "lalpha")

  # Vérification des colonnes
  expect_true(all(c("vario", "chi") %in% colnames(result)))

  # Vérification des dimensions
  expect_equal(nrow(result), nrow(df_lags))

  # Vérification des valeurs de chi (bornées entre 1e-8 et 1)
  expect_true(all(result$chi >= 0 & result$chi <= 1))
})

test_that("theoretical_chi computes correct variogram values", {

  df_lags <- data.frame(
    s1 = c(1),
    s2 = c(2),
    tau = c(1),
    s1x = c(0),
    s1y = c(0),
    s2x = c(1),
    s2y = c(1),
    hnorm = c(1.41)
  )

  # No advection, no wind
  params <- c(0.5, 0.3, 1.2, 0.8, 0, 0)
  result <- theoretical_chi(params, df_lags, latlon = FALSE,
                            distance = "euclidean")
  # Verif variogram value
  expected_vario <- 2 * params[1] * abs(df_lags$hnorm)^params[3] +
                    2 * params[2] * abs(df_lags$tau)^params[4]
  expect_equal(result$vario, expected_vario, tolerance = 1e-2)

  # Verif chi value
  chi_value <- 2 * (1 - pnorm(sqrt(0.5 * expected_vario)))
  chi_value <- max(min(chi_value, 1), 1e-8)  # Bornage
  expect_equal(result$chi, chi_value, tolerance = 1e-2)
})



test_that("theoretical_chi", {
  tau_vect <- 0:10
  data_rain <- matrix(runif(100 * 25, 0, 100), nrow = 100, ncol = 25)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(5) # Example sites_coords
  true_param <- c(0.4, 0.2, 1.5, 1, 0, 0)
  df_lags <- get_lag_vectors(sites_coords, tau_vect = -10:10)

  chi_theorical <- theoretical_chi(true_param, df_lags, distance = "euclidean")

  # for one hnorm and one tau
  hnorm <- 1
  tau <- 3
  semivar <- true_param[1]*hnorm^true_param[3] + true_param[2]*tau^true_param[4]
  chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

  chi_h_t <- chi_theorical$chi[chi_theorical$hnorm == hnorm &
                                    chi_theorical$tau == tau]

  expect_equal(unique(chi_h_t), chi_h_t_verif)

  # for r-pareto, conditionning on s0, t0
  s0 <- c(1, 1)
  t0 <- 1
  adv <- c(0.5, 0.3)
  true_param <- c(true_param, adv)

  df_lags <- get_conditional_lag_vectors(sites_coords, tau_vect = 0:10,
                                         s0 = s0, t0 = t0)

  chi_theorical <- theoretical_chi(true_param, df_lags, distance = "euclidean")

  # for one hnorm and one tau
  hnorm <- chi_theorical$hnorm[10]
  tau <- chi_theorical$tau[10]
  semivar <- true_param[1]*hnorm^true_param[3] + true_param[2]*tau^true_param[4]
  chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

  chi_h_t <- chi_theorical$chi[chi_theorical$hnorm == hnorm &
                                    chi_theorical$tau == tau]

  expect_equal(unique(chi_h_t), chi_h_t_verif)
})


test_that("theoretical_chi computes correct chi values with complex conditions", {
  set.seed(42)  # For reproducibility
  # Create random and irregular data
  tau_vect <- 0:10
  data_rain <- matrix(runif(100 * 30, 0, 100), nrow = 100, ncol = 30)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(30)  # More complex grid

  # Variogram parameters
  true_param <- c(0.5, 0.25, 1.7, 1.2, 0, 0)  # beta1, beta2, alpha1, alpha2,
                                              # adv_x, adv_y
  df_lags <- get_lag_vectors(sites_coords, tau_vect = -10:10)

  # Calculate theoretical chi without advection
  chi_theorical <- theoretical_chi(true_param, df_lags, distance = "euclidean")

  # Comparison for multiple (hnorm, tau) values
  for (i in sample(1:nrow(df_lags), 5)) {  # Test on 5 random values
    hnorm <- df_lags$hnorm[i]
    tau <- df_lags$tau[i]
    semivar <- true_param[1] * hnorm^true_param[3] + true_param[2] *
                                                        abs(tau)^true_param[4]
    chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

    chi_h_t <- chi_theorical$chi[chi_theorical$hnorm == hnorm &
                                    chi_theorical$tau == tau]

    expect_equal(unique(chi_h_t), chi_h_t_verif, tolerance = 1e-6)
  }

  # **Case with conditioning on s0, t0 and advection**
  s0 <- c(2, 3)  # Randomly chosen reference site
  t0 <- 2
  adv <- c(0.4, 0.2)  # Advection in x and y
  true_param <- c(0.5, 0.25, 1.7, 1.2, adv)

  df_lags <- get_conditional_lag_vectors(sites_coords, tau_vect = 0:10,
                                          s0 = s0, t0 = t0)
})


test_that("theoretical_chi correctly accounts for wind effect", {
  set.seed(42)  # Pour reproductibilité
  
  # Création d'une grille plus complexe
  tau_vect <- 0:10
  data_rain <- matrix(runif(100 * 30, 0, 100), nrow = 100, ncol = 30)
  data_rain <- as.data.frame(data_rain)
  sites_coords <- generate_grid_coords(30)  # Grille irrégulière

  # Paramètres du variogramme avec wind_vect
  true_param <- c(0.5, 0.25, 1.7, 1.2, 0.8, 0.3)  # beta1, beta2, alpha1, alpha2, eta1, eta2
  df_lags <- get_lag_vectors(sites_coords, tau_vect = -10:10)

  # Cas sans vent
  chi_theorical_no_wind <- theoretical_chi(true_param, df_lags,
                                                  distance = "euclidean")

  # Création d'un vecteur de vent aléatoire par site
  wind_vect <- c(-2, 3)

  # Calcul du chi théorique avec vent
  chi_theorical_wind <- theoretical_chi(true_param, df_lags, distance = "euclidean")

  # Test sur plusieurs valeurs (hnorm, tau)
  for (i in sample(1:nrow(df_lags), 5)) {
    hnorm <- df_lags$hnorm[i]
    tau <- df_lags$tau[i]

    # Effet du vent (advection dépendant du vent)
    wind_kmh <- wind_vect * 3.6  # Conversion en km/h
    adv_x <- abs(wind_kmh)^true_param[5] * sign(wind_kmh) * true_param[6]

    # Recalcul de la distance avec l'effet du vent
    s2xv <- df_lags$s2x[i] - adv_x * tau  # Déplacement du site s2 selon le vent
    hnormV <- sqrt((s2xv - df_lags$s1x[i])^2 + (df_lags$s2y[i] -
                                                            df_lags$s1y[i])^2)

    semivar <- true_param[1] * hnormV^true_param[3] + true_param[2] *
                                                          abs(tau)^true_param[4]
    chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

    chi_h_t_wind <- chi_theorical_wind$chi[chi_theorical_wind$hnorm == hnorm &
                                           chi_theorical_wind$tau == tau]
    expect_equal(max(unique(chi_h_t_wind), 1e-8), max(chi_h_t_verif, 1e-8),
                                                            tolerance = 1e-3)
  }
})



test_that("theoretical_chi correctly accounts for wind effect in lat/lon", {
  set.seed(42)  # Pour reproductibilité
  
  # Création d'une grille avec des latitudes et longitudes réalistes
  sites_coords <- generate_realistic_latlon_grid(30)  # Grille irrégulière avec lat/lon
  
  
  df_lags <- get_lag_vectors(sites_coords, tau_vect = -10:10, latlon = TRUE)
  df_lags <- df_lags %>%
    dplyr::filter(tau >= -10 & tau <= 10)  # On garde tau dans [-10,10]

  # Paramètres du variogramme avec le vent
  eta1 <- 0.8
  eta2 <- 0.3
  true_param <- c(0.5, 0.25, 1.7, 1.2)  # beta1, beta2, alpha1, alpha2

  # Création d'un vecteur de vent aléatoire en m/s par site
  wind_vect <- c(-2, 3)
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  true_param <- c(true_param, adv)  # Ajout de l'advection aux paramètres

  # Calcul du chi théorique avec vent
  chi_theorical_wind <- theoretical_chi(true_param, df_lags, latlon = TRUE,
                                        distance = "euclidean")

  # Test sur plusieurs valeurs (hnorm, tau) avec coordonnées géographiques
  for (i in sample(1:nrow(df_lags), 5)) {
    hnorm <- df_lags$hnorm[i]
    tau <- df_lags$tau[i]

    # Conversion du vent en advection avec eta1 et eta2
    wind_kmh <- wind_vect[i] * 3.6  # Conversion en km/h
    adv_x <- abs(wind_kmh)^true_param[5] * sign(wind_kmh) * true_param[6]  

    # Conversion des coordonnées (lat/lon → km)
    lat_convert <- 111.32  # km par degré de latitude
    lon_convert <- lat_convert * cos(pi * df_lags$s2y[i] / 180)  # Correction longitude

    # Application de l'effet du vent en tenant compte de la conversion
    s2xv <- df_lags$s2x[i] * lon_convert - adv_x * tau
    s2yv <- df_lags$s2y[i] * lat_convert

    # Recalcul de la distance avec l'effet du vent
    hnormV <- sqrt((s2xv - df_lags$s1x[i] * lon_convert)^2 + (s2yv - df_lags$s1y[i] * lat_convert)^2)

    # Calcul de la semi-variogramme
    semivar <- true_param[1] * hnormV^true_param[3] + true_param[2] * abs(tau)^true_param[4]
    chi_h_t_verif <- 2 * (1 - pnorm(sqrt(semivar)))

    # Vérification de chi
    chi_h_t_wind <- chi_theorical_wind$chi[chi_theorical_wind$hnormV == hnormV &
                                           chi_theorical_wind$tau == tau]

    expect_equal(unique(chi_h_t_wind), chi_h_t_verif, tolerance = 1e-6)
  }
})


test_that("Latitude/Longitude conversion to meters is correct", {
  # Known coordinates (near the equator)
  lat1 <- 0
  lon1 <- 0
  lat2 <- 0
  lon2 <- 1  # 1 degree to the east

  # Expected distance (in meters)
  expected_distance_km <- 111.195  # 1 degree longitude at equator
  expected_distance_m <- expected_distance_km * 1000  # Convert to meters

  sites_coords <- data.frame(
    site = 1:2,
    Longitude = c(lon1, lon2),
    Latitude = c(lat1, lat2)
  )

  df_lags <- get_lag_vectors(sites_coords, tau_vect = 0, latlon = TRUE)
  df_lags$hnorm <- df_lags$hnorm # Convert to km
  expect_true(round(df_lags$hnorm[1], 3) == expected_distance_km)
  expect_true(round(df_lags$hnorm[2], 2) == 0)
  # Assuming your function for conversion is like this:
  # Convert degrees to kilometers (latlon = TRUE)
  lat_convert <- 111.32  # km per degree latitude (this is constant at the equator)
  lon_convert <- lat_convert * cos(pi * df_lags$s2y / 180)  # longitude conversion factor

  params <- c(0.5, 0.25, 1.7, 1.2, 0.8, 0.3)  # Example variogram parameters
  chi_df <- theoretical_chi(params, df_lags, latlon = TRUE, distance = "lalpha")
  # Check that the calculated distance is close to the expected distance
  expect_equal(chi_df$hlag[1], expected_distance_km, tolerance = 1e-1)  # Allow small tolerance


  # Known coordinates (near the equator)
  lat1 <- 43.62505
  lon1 <- 3.862038 # Montpellier
  lat2 <- 43.683333
  lon2 <- 4.133333  # Lunel

  sites_coords <- data.frame(
    site = 1:2,
    Longitude = c(lon1, lon2),
    Latitude = c(lat1, lat2)
  )

  # Expected distance (in meters)
  expected_distance_km <- 22.77  # See on internet
  expected_distance_m <- expected_distance_km * 1000  # Convert to meters


  # Create a dataframe with these coordinates
  df_lags <- get_lag_vectors(sites_coords, tau_vect = 0:10, latlon = TRUE)
  df_lags$hnorm <- df_lags$hnorm # Convert to km
  expect_true(round(df_lags$hnorm[1], 2) == expected_distance_km)

  eta1 <- 0.8
  eta2 <- 0.3
  params <- c(0.5, 0.25, 1.7, 1.2) # Example variogram parameters
  wind_vect <- c(10, 5)  # Example wind vector (vx = 10 m/s, vy = 5 m/s)
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  params <- c(params, adv)  # Add advection to the parameters
  
  chi_df <- theoretical_chi(params, df_lags,
                            latlon = TRUE, distance = "euclidean")

  # hnormV change with the wind for each tau
  expect_false(chi_df$hnormV[1] == chi_df$hnormV[2]) # tau = 0 and tau = 1 for
                                                     # the same site pair

})

test_that("theoretical_chi computes correct values with real lat/lon coordinates, wind_vect", {
  # Sample parameters
  eta1 <- 0.5
  eta2 <- 0.5
  params <- c(1, 1, 1, 1)

  # Real-world lat/lon coordinates for Paris and London
  df_lags_PL <- data.frame(
    s1 = c(1), s2 = c(2), tau = c(0),
    s1x = c(2.3522), s1y = c(48.8566),  # Paris (longitude, latitude)
    s2x = c(-0.1278), s2y = c(51.5074), # London (longitude, latitude)
    hnorm = c(343.51)  # Approx distance in km
  )

  # Wind vector
  wind_vect <- c(0, 0)  # No wind effect for a simple test case
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  params <- c(params, adv)  # Add advection to the parameters

  # Compute theoretical chi with latlon = TRUE
  result <- theoretical_chi(params, df_lags_PL, latlon = TRUE,
                          distance = "euclidean")

  # Expected values assuming haversine distance with 0° at North
  expected_theta <- atan2(51.5074 - 48.8566, -0.1278 - 2.3522)
  expected_x_polar <- df_lags_PL$hnorm * cos(expected_theta)
  expected_y_polar <- df_lags_PL$hnorm * sin(expected_theta)

  # Check theta, x_polar, and y_polar values
  expect_equal(result$theta, expected_theta, tolerance = 1e-6)
  expect_equal(round(result$hnormV, 0), round(df_lags_PL$hnorm, 0))
  expect_equal(result$x_polar, expected_x_polar, tolerance = 1e-1)
  expect_equal(result$y_polar, expected_y_polar, tolerance = 1e-1)


  df_lags_LP <- data.frame(
    s1 = c(1), s2 = c(2), tau = c(0),
    s1x = c(2.3522), s1y = c(48.8566),  # Paris (longitude, latitude)
    s2x = c(-0.1278), s2y = c(51.5074), # London (longitude, latitude)
    hnorm = c(343.51)  # Approx distance in km
  )

  
  # Compute theoretical chi with latlon = TRUE
  result <- theoretical_chi(params, df_lags_LP, latlon = TRUE,
                           distance = "lalpha")

  expected_theta <- atan2(51.5074 - 48.8566, -0.1278 - 2.3522)
  # Expected values assuming haversine distance with 0° at North
  expected_theta_meteo <- (450 - atan2(51.5074 - 48.8566, -0.1278 - 2.3522) *
                                                              180 / pi) %% 360
  expected_x_polar <- df_lags_LP$hnorm * cos(expected_theta)
  expected_y_polar <- df_lags_LP$hnorm * sin(expected_theta)

  expect_equal(result$theta, expected_theta, tolerance = 1e-6)
  expect_equal(round(result$hnormV, 0), round(df_lags_LP$hnorm, 0))
})




test_that("theoretical_chi computes correct values with different etas", {
  # Sample parameters
  eta1 <- 0.5
  eta2 <- 0.5
  params <- c(1, 1, 1, 1)

  # Real-world lat/lon coordinates for Paris and London
  df_lags <- data.frame(
    s1 = c(1), s2 = c(2), tau = c(0),
    s1x = c(2.3522), s1y = c(48.8566),  # Paris (longitude, latitude)
    s2x = c(2.3523), s2y = c(48.8566), # Close to Paris (longitude, latitude)
    hnorm = c(0.1)  # Approx distance in km
  )

  # Create new dataframe with tau = 0 and tau = 1
  df_lags_extended <- df_lags[rep(1, 2), ]
  df_lags_extended$tau <- c(1,2)

  # Wind vector
  wind_vect <- c(-9.8, 8)  # Wind effect for a simple test case
  adv <- eta1 * sign(wind_vect) * abs(wind_vect)^eta2
  params <- c(params, adv)  # Add advection to the parameters

  params_noadv <- c(1, 1, 1, 1, 0, 0)
  result_noadv <- theoretical_chi(params_noadv, df_lags_extended, latlon = TRUE,
                          distance = "lalpha")

  # Compute theoretical chi with latlon = TRUE
  result_1 <- theoretical_chi(params, df_lags_extended, latlon = TRUE,
                           distance = "lalpha")

  params_2 <- c(1, 1, 1, 1, -1, -1)
  result_2 <- theoretical_chi(params_2, df_lags_extended, latlon = TRUE,
                          distance = "lalpha")

  # Check that the chi values are different for different etas
  expect_true(any(result_1$chi != result_2$chi))

})
