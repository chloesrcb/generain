# --- Paramètres utilisateur / hyperparam ---
max_iter <- 10
tol_eta_rel <- 1e-2           # tolérance relative sur eta pour converger
init_eta <- c(1, 1)           # démarrage
hmax_quantile_choice <- 0.9   # quantile pour proposer hmax à chaque itération (0.75/0.9/0.95)
verbose <- TRUE

# helper: calcule adv (en km/s si wind en m/s et on veut hnormV en km)
make_adv_mat <- function(wind_df, eta1, eta2) {
  adv_x <- (abs(wind_df$vx)^eta2) * sign(wind_df$vx) * eta1
  adv_y <- (abs(wind_df$vy)^eta2) * sign(wind_df$vy) * eta1
  cbind(adv_x, adv_y)
}

# helper: compute hnormV vector for df_lags (row-wise adv)
compute_hnormV <- function(df_lags, adv_mat) {
  # attente: adv_mat a autant de lignes que df_lags (ou est broadcastable)
  tau <- df_lags$tau
  s1x <- df_lags$s1x; s1y <- df_lags$s1y
  s2x <- df_lags$s2x; s2y <- df_lags$s2y
  # if adv_mat is per-episode, make sure to align it to df_lags rows externally
  s2xv <- s2x - adv_mat[,1] * tau
  s2yv <- s2y - adv_mat[,2] * tau
  sqrt((s2xv - s1x)^2 + (s2yv - s1y)^2)
}

# main iterative wrapper
iterative_hmax_eta <- function(df_lags, wind_df, params_init4,
                               list_episodes, list_excesses, list_lags,
                               neg_ll_composite,
                               max_iter = 8, tol = 1e-2,
                               quantile_hmax = 0.9, clip_wind = NULL,
                               verbose = TRUE) {

  eta_old <- init_eta
  params_best <- c(params_init4, eta_old) # initial full params
  for (it in 1:max_iter) {
    if (verbose) cat("Iter", it, " - eta_old =", paste(round(eta_old,4), collapse=", "), "\n")

    # optionnel: clip winds pour diagnostic (évite influence extrême)
    wind_tmp <- wind_df
    if (!is.null(clip_wind)) {
      wind_tmp$vx <- pmax(pmin(wind_tmp$vx, clip_wind), -clip_wind)
      wind_tmp$vy <- pmax(pmin(wind_tmp$vy, clip_wind), -clip_wind)
    }

    # 1) compute adv with current eta and hnormV
    adv_mat <- make_adv_mat(wind_tmp, eta_old[1], eta_old[2])
    # Ensure adv_mat has same nrows as df_lags (align if adv is per-episode)
    if (nrow(adv_mat) != nrow(df_lags)) {
      # try to recycle per-episode adv to rows (user must ensure mapping)
      if (nrow(adv_mat) == length(list_episodes)) {
        # assume one adv per episode, expand to df_lags via an 'episode' column
        if (!"episode" %in% names(df_lags)) stop("df_lags must have 'episode' column if adv per-episode")
        adv_mat_rows <- adv_mat[df_lags$episode, , drop = FALSE]
      } else {
        stop("adv_mat rows mismatch df_lags rows; adapt alignment")
      }
    } else adv_mat_rows <- adv_mat

    hnormV_vec <- compute_hnormV(df_lags, adv_mat_rows)

    # propose hmax based on quantile
    hmax_prop <- quantile(hnormV_vec, probs = quantile_hmax, na.rm = TRUE)
    if (verbose) cat(" Proposed hmax (q=", quantile_hmax, ") =", round(hmax_prop,4), "\n")

    # 2) estimate parameters with current hmax (pass fixed_eta NA to estimate)
    # Build initial param vector for optimizer (beta1,beta2,alpha1,alpha2, eta1,eta2)
    params_start <- c(params_init4, eta_old)
    # call optimizer (adapt to your optimizer; here oneliner with optim as example)
    opt_fun <- function(par) {
      # pass hmax to neg_ll_composite, and wind_df so adv computed inside
      neg_ll_composite(params = par, list_episodes = list_episodes,
                       list_excesses = list_excesses, list_lags = list_lags,
                       wind_df = wind_df, hmax = hmax_prop, latlon = TRUE,
                       distance = "euclidean", threshold = FALSE, rpar = TRUE)
    }
    # bounds and transforms recommended (not shown). Use optim as example:
    res_opt <- optim(par = params_start, fn = opt_fun, method = "L-BFGS-B",
                     lower = c(1e-6, 1e-6, 0.01, 0.01, 1e-6, 1e-6),
                     upper = c(Inf, Inf, 2, 2, 1e3, 1e3),
                     control = list(maxit = 500))
    params_est <- res_opt$par
    eta_new <- params_est[5:6]
    if (verbose) cat(" Estimated eta_new =", paste(round(eta_new,4), collapse=", "), "\n")

    # 3) convergence check
    rel_change <- abs(eta_new - eta_old) / pmax(abs(eta_old), 1e-8)
    if (verbose) cat(" rel change eta:", paste(round(rel_change,4), collapse=", "), "\n")
    params_best <- params_est

    if (all(rel_change < tol)) {
      if (verbose) cat("Converged (tol reached).\n")
      break
    }

    # update for next iter
    eta_old <- eta_new
    # update init 4 params to speed convergence (optional)
    params_init4 <- params_est[1:4]
  }

  list(params = params_best, iters = it, converged = all(rel_change < tol), hmax_final = hmax_prop, opt_res = res_opt)
}

# --- Example call (adapte les objets ci-dessous) ---
# params_init4 = c(beta1, beta2, alpha1, alpha2) initial guess
params_init4 <- c(0.2, 1.0, 0.65, 0.73)

df_lags_list <- list()

for (i in seq_along(list_lags)) {
  tmp <- list_lags[[i]]
  tmp$episode <- i   # ajoute un identifiant épisode
  df_lags_list[[i]] <- tmp
}

new_df_lags <- do.call(rbind, df_lags_list)

res_iter <- iterative_hmax_eta(
  df_lags = new_df_lags,
  wind_df = wind_opt,
  params_init4 = params_init4,
  list_episodes = list_episodes,
  list_excesses = list_excesses,
  list_lags = list_lags,
  neg_ll_composite = neg_ll_composite,
  max_iter = 8,
  tol = 1e-2,
  quantile_hmax = 0.9,
  clip_wind = 10,            # optional clipping at 10 m/s for diagnostics
  verbose = TRUE
)

print(res_iter)
