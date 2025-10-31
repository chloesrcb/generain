chi <- theoretical_chi(params, df_lags, latlon, distance)
# Choisir la colonne de distance pour troncature (advectée si dispo)
h_col <- if ("hnormV" %in% names(chi)) "hnormV" else "hnorm"

# (A) Couper aux h où chi est devenu négligeable (ex: 0.05 ou 0.1)
chi_cut <- 0.05
chi <- as.data.frame(chi)
summary(chi$chi)
min(chi$chi)
chi_no1 <- chi$chi[chi$chi < 0.99]
max(chi_no1)
# Option 1 : seuil sur chi uniquement
chi_keep_A <- chi[chi$chi >= chi_cut, ]
length(chi$chi)
length(chi_keep_A$chi)
# Option 2 : seuil mixte chi & distance (souvent plus stable)
hmax_guess <- with(chi, tapply(get(h_col), tau, function(v) quantile(v, 0.9, na.rm=TRUE)))
# sinon fixe: hmax_guess <- 60   # ex. 60 km
chi_keep_B <- chi[chi$chi >= chi_cut & chi[[h_col]] <= hmax_guess[as.character(chi$tau)], ]

# médiane(chi) par tau et par classes de distance
aggregate(chi ~ tau, data = chi, median)
chi_df <- as.data.frame(chi)
chi_df$h_bin <- cut(chi_df[[h_col]], breaks = 10)
aggregate(chi ~ h_bin, data = chi_df, median)
table(close_raw = chi$hnorm <= 40, close_adv = chi$hnormV <= 40)
hmax <- 30  # km
tau_max <- 8
chi_keep <- chi[chi[[h_col]] <= hmax & chi$tau <= tau_max, ]
chi$h_bin <- cut(chi_df[[h_col]], breaks = seq(0, max(chi[[h_col]]), by = 5))
aggregate(chi ~ h_bin, data = chi_df, median)
library(ggplot2)
ggplot(chi_df, aes(x = hnorm, y = chi, color = factor(tau))) +
  geom_point(alpha = 0.5) +
  geom_smooth(se = FALSE) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = "Distance effective (hnormV, km)", y = "Chi",
       color = "Tau (heures)",
       title = "Décroissance spatiale et temporelle de Chi")