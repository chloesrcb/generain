# Example grid dimensions
nx <- 5  # number of points in x dimension
ny <- 5  # number of points in y dimension
nt <- 3  # number of points in t dimension

# Create example grid coordinates
x <- seq(1, 10, length.out = nx)
y <- seq(1, 10, length.out = ny)
t <- seq(0.1, 10, length.out = nt)

# Create spatio-temporal grid
grid <- expand.grid(x = x, y = y, t = t)

# Extract coordinates
coords <- grid[, c("x", "y", "t")]

# Calculate Euclidean distances
distances <- dist(coords)

model <- RMfbm(alpha=1)
x <- seq(0, 10, 0.02)
plot(model)
plot(RFsimulate(model, x=x))

modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) + 
                 RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 2)
x <- seq(0, 10, 0.02)
t <- seq(0, 10, 0.02)

plot(RFsimulate(modelBuhlCklu, x=x, T=8, seed=1))


modelBuhlCklu <- RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 1) + 
                 RandomFields::RMfbm(alpha = alpha1, var = beta1, proj = 2) +
                 RandomFields::RMfbm(alpha = alpha2, var = beta2, proj = 3)
x <- seq(0, 10, 0.02)
y <- seq(0, 10, 0.02)
t <- seq(0, 10, 0.02)
plot(RFsimulate(modelBuhlCklu, x=x, y=1:10, T=8, seed=1))


W <- RandomFields::RFsimulate(model, x = x, y = y, T = t,
                             distances = distances, dim=2)


# Example grid of coordinates
x <- seq(0, 10, by = 1)
y <- seq(0, 10, by = 1)
t <- seq(0, 1, by = 0.1)

# Construct grid
grid <- expand.grid(x = x, y = y, t = t)

# Compute distances
dx <- outer(x, x, "-")
dy <- outer(y, y, "-")
dt <- outer(t, t, "-")
norm_h <- sqrt(dx^2 + dy^2)  # Spatial distances
distances <- dist(coords)
W <- RFsimulate(model = modelBuhlCklu, distances = distances, dim=3)
