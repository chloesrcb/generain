library(generain)
library(reshape2)
library(ggplot2)
setwd("./script")
# load libraries
source("load_libraries.R")

library(animation)
library(ggplot2)
library(plotly)

# spatial and temporal structures
ngrid <- 10
spa <- 1:ngrid
nsites <- ngrid^2 # if the grid is squared
temp <- 1:100

# beta1, beta2, alpha1, alpha2
param <- c(0.8, 0.4, 1.5, 1) # true parameters for the variogram
beta1 <- param[1] / 2
beta2 <- param[2] / 2
# beta3 <- param[2] / 2
alpha1 <- param[3]
alpha2 <- param[4]


file_path <- paste0("../data/simulations_BR/sim_", nsites,
          "s_",  length(temp), "t/br_", nsites, "s_",  length(temp), "t_1.csv")

simulation_data <- read.csv(file_path)
ngrid <- sqrt(ncol(simulation_data))  # Number of grid points in each dimension

simulation_data$Time <- rownames(simulation_data) # Add a time column
simulation_data_long <- melt(simulation_data) # Convert to long format
simulation_data_long$Time <- as.numeric(simulation_data_long$Time)
# Create a dataframe to represent grid points
grid <- expand.grid(x = 1:ngrid, y = 1:ngrid)

plots <- list()
# for each time step
for (i in unique(simulation_data_long$Time)) {
  # Add the simulated values to the grid dataframe
  grid$value <- simulation_data_long$value[simulation_data_long$Time == i]

  # Plot
  p <-  ggplot(data = grid, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_gradient(low = "#70a7ae", high = "#9d503d",
                        name = "Rainfall in mm",
                        limits = c(min(simulation_data_long$value),
                        max(simulation_data_long$value))
                       ) +
    labs(title = paste0("t =", i, " | Betas: ", beta1, ", ", beta2,
                        " | Alphas: ",
                        alpha1, ", ", alpha2, " | Advection: ", adv[1],
                        ", ", adv[2])) +
    theme_minimal() +
    theme(plot.background = element_rect(fill = "#F9F8F6", color = "#F9F8F6"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())

  plots[[i]] <- p
}

# Save the plots as a gif
ani.options(interval = 0.5) # time between frames
saveGIF({
  for (i in 1:ntimes) {
    print(plots[[i]])
  }
}, movie.name = paste0("/user/cserreco/home/Documents/These/generain/images",
                       "/simu_gif/sim_br/br_", ngrid^2, "s_",
                                length(temp), ".gif"),
    ani.width = 700, ani.height = 600, ani.units = "px", ani.type = "cairo")


data <- data.frame(
    x = grid$x,
    y = grid$y,
    t = simulation_data_long$Time,
    value = as.vector(simulation_data_long$value[simulation_data_long$Time == t])
)



spa <- 1:10
selected_time <- 1
rain_data <- simulation_data_long[simulation_data_long$Time == selected_time,]
grid_data <- expand.grid(x = spa, y = spa)
rain_matrix <- matrix(rain_data$value, nrow = length(spa), ncol = length(spa))
# grid_data$z <- as.vector(rain_data$value)

fig <- plot_ly(
  x = ~spa,
  y = ~spa,
  z = ~rain_matrix,
  type = 'surface',
  colorscale = 'Viridis',
  zmin = 0, 
  zmax = max(rain_matrix), 
  colorbar = list(title = "Rainfall amount (mm)")
)

fig <- fig %>% layout(
  scene = list(
    xaxis = list(title = "x"),
    yaxis = list(title = "y"),
    zaxis = list(title = "Rainfall (mm)")
  )
)

fig



lon <- 1:10
lat <- 1:10

axx <- list(

  gridcolor='rgb(255, 255, 255)',

  zerolinecolor='rgb(255, 255, 255)',

  showbackground=TRUE,

  backgroundcolor='rgb(230, 230,230)'

)


# individual plots
spa <- 1:10
selected_time <- 1
rain_data <- simulation_data_long[simulation_data_long$Time == selected_time,]
grid_data <- expand.grid(x = spa, y = spa)
rain_matrix1 <- matrix(rain_data$value, nrow = length(spa), ncol = length(spa))

fig1 <- plot_ly(
  x = ~spa,
  y = ~spa,
  z = ~rain_matrix1,
  scene='scene1'
)

fig1 <- fig1 %>% add_surface(showscale=FALSE)


selected_time <- 2
rain_data <- simulation_data_long[simulation_data_long$Time == selected_time,]
grid_data <- expand.grid(x = spa, y = spa)
rain_matrix2 <- matrix(rain_data$value, nrow = length(spa), ncol = length(spa))

fig2 <- plot_ly(
  x = ~spa,
  y = ~spa,
  z = ~rain_matrix2,
  scene='scene2'
)

fig2 <- fig2 %>% add_surface(showscale=FALSE)


selected_time <- 3
rain_data <- simulation_data_long[simulation_data_long$Time == selected_time,]
grid_data <- expand.grid(x = spa, y = spa)
rain_matrix3 <- matrix(rain_data$value, nrow = length(spa), ncol = length(spa))

fig3 <- plot_ly(
  x = ~spa,
  y = ~spa,
  z = ~rain_matrix3,
  scene='scene3'
)

fig3 <- fig3 %>% add_surface(showscale=FALSE)


selected_time <- 4
rain_data <- simulation_data_long[simulation_data_long$Time == selected_time,]
grid_data <- expand.grid(x = spa, y = spa)
rain_matrix4 <- matrix(rain_data$value, nrow = length(spa), ncol = length(spa))

fig4 <- plot_ly(
  x = ~spa,
  y = ~spa,
  z = ~rain_matrix4,
  scene='scene4'
)

fig4 <- fig4 %>% add_surface(showscale=FALSE)


# subplot and define scene

fig <- subplot(fig1, fig2, fig3, fig4) 

fig <- fig %>% layout(title = "3D Subplots",

         scene = list(domain=list(x=c(0,0.5),y=c(0.5,1)),

                      xaxis=axx, yaxis=axx, zaxis=axx,

                      aspectmode='cube'),

         scene2 = list(domain=list(x=c(0.5,1),y=c(0.5,1)),

                       xaxis=axx, yaxis=axx, zaxis=axx,

                       aspectmode='cube'),

         scene3 = list(domain=list(x=c(0,0.5),y=c(0,0.5)),

                       xaxis=axx, yaxis=axx, zaxis=axx,

                       aspectmode='cube'),

         scene4 = list(domain=list(x=c(0.5,1),y=c(0,0.5)),

                       xaxis=axx, yaxis=axx, zaxis=axx,

                       aspectmode='cube'))


fig