
# compare sim_rpareto vs sim_rpareto_coords
lx <- 5
ly <- 5
x0 <- 1
y0 <- 1
df_coords <- expand.grid(
  Longitude = 1:lx,
  Latitude  = 1:ly
)
df_coords$site <- paste0("S", seq_len(nrow(df_coords)))
s0_site <- with(df_coords, site[Longitude==x0 & Latitude==y0][1])
s0_coords <- with(df_coords, c(Longitude[site == s0_site][1], Latitude[site == s0_site][1]))
s0_index <- which(df_coords$site == s0_site)


set.seed(123)
Z_rpar <- sim_rpareto(
  beta1 = 1, beta2 = 1, alpha1 = 1, alpha2 = 1, 
  x= 1:lx, y = 1:ly,
  t = 0:11, # every 5 minutes for 1 hour
    s0 = s0_coords, t0 = 0, nres = 1, seed = 123
)
dim(Z_rpar$Z) # 5 x 5 x 12 x 1 x 1

Z_sim_rpar <- Z_rpar$Z[, , , 1, 1] # 5 x 5 x 12

Z_rpar_coords <- sim_rpareto_coords(
  beta1 = 1, beta2 = 1, alpha1 = 1, alpha2 = 1, 
  coords = df_coords[, c("Longitude", "Latitude")], 
  times = 0:11, # every 5 minutes for 1 hour
    s0_index = s0_index, t0_index = 1, seed = 123
)

dim(Z_rpar_coords$Z) # 25 x 12
Z_sim_rpar_coords <- Z_rpar_coords$Z

Z_rpar_sites <- sim_rpareto_sites(
  beta1 = 1, beta2 = 1, alpha1 = 1, alpha2 = 1, 
  coords = df_coords, 
  t = 0:11, # every 5 minutes for 1 hour
    s0 = s0_site, t0 = 0, seed = 123
)

dim(Z_rpar_sites$Z) # 25 x 12
Z_sim_rpar_sites <- Z_rpar_sites$Z

# compare
Z_sim_rpar[1, 1, ] # 12 values for site 1
Z_sim_rpar_coords[1, ] # 12 values for site 1
Z_sim_rpar_sites[1, ] # 12 values for site 1

# check cumul by site
apply(Z_sim_rpar, c(1, 2), sum) # 5 x 5
apply(Z_sim_rpar_coords, 1, sum) # 25
apply(Z_sim_rpar_sites, 1, sum) # 25

