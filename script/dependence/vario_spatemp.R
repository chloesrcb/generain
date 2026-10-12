################################################################################
## Bootstrap IC pforour gamma(h,tau)-gamma(0,tau)
################################################################################

library(dplyr)
library(purrr)

dates <- as.POSIXct(rownames(rain_new))
month_id <- format(dates,"%Y-%m")
months <- unique(month_id)

B <- 10 

boot_sep <- vector("list", B)

for(b in 1:B){

  cat("Bootstrap", b, "/", B, "\n")
  sampled_months <- sample(months,
                           length(months),
                           replace = TRUE)

  ind <- unlist(lapply(sampled_months,
                       function(m) which(month_id==m)))

  rain_boot <- rain_new[ind,]

  ## chi
  chi_boot <- calc_chi_spatemp_with_zero(
    rain_boot,
    dist_mat,
    tau_max = 10,
    q = q,
    remove_zeros = FALSE
  )

  chi_boot <- chi_boot %>%
    mutate(
      dist_class = cut(distance,
                       breaks = radius_posits,
                       include.lowest = TRUE)
    ) %>%
    mutate(
      dist_class = as.character(dist_class),
      dist_class = ifelse(distance==0,"0",dist_class)
    )

  ## mean by dist_class and tau
  chi_summary_boot <- chi_boot %>%
    group_by(dist_class,tau) %>%
    summarise(
      mean_chi = mean(chi,na.rm=TRUE),
      .groups="drop"
    ) %>%
    left_join(df_prev_h,by="dist_class")

  ## variogram
  vario_boot <- chi_summary_boot %>%
    mutate(vario = vario_spatemp(mean_chi))

  ## reference h=0
  ref0 <- vario_boot %>%
    filter(h==0) %>%
    select(tau,vario0=vario)

  ## gamma(h,tau)-gamma(0,tau)
  boot_sep[[b]] <- vario_boot %>%
    left_join(ref0,by="tau") %>%
    mutate(diff = vario-vario0,
           boot=b)
}

boot_sep <- bind_rows(boot_sep)

IC <- boot_sep %>%
  filter(h>0) %>%
  group_by(h,tau) %>%
  summarise(
    lower = quantile(diff,0.025,na.rm=TRUE),
    upper = quantile(diff,0.975,na.rm=TRUE),
    .groups="drop"
  )


df_plot <- df_separabilite_zero %>%
  filter(h>0) %>%
  left_join(IC,by=c("h","tau"))

df_plot <- df_plot %>%
  mutate(
    tau_min = factor(5 * tau,
                     levels = sort(unique(5 * tau)))
  )

ggplot(df_plot,
       aes(x = h,
           y = diff_vario_0,
           colour = tau_min,
           fill = tau_min)) +
  geom_ribbon(aes(ymin=lower,
                  ymax=upper),
              alpha=0.20,
              colour=NA) +
  geom_smooth(method="loess",
              se=FALSE,
              span=3,
              linewidth=1) +
  labs(x = "Distance h (m)",
    y = expression(hat(gamma)(bold(h),tau)-hat(gamma)(0,tau)),
    colour = "Time lag (min)",
    fill = "Time lag (min)"
    ) +
  btf_theme
