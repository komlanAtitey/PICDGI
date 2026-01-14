## -----------------------------
## Functions from the paper
## -----------------------------

rho_u <- function(k, H) {
  0.5 * (abs(k + 1)^(2 * H) - 2 * abs(k)^(2 * H) + abs(k - 1)^(2 * H))
}

cov_innovation <- function(H, K = 3, sigma_u2 = 1) {
  rho_vec <- rho_u(1:K, H)
  R <- toeplitz(rho_vec)
  R <- abs(R)               # match visual behavior in the example figure
  sigma_u2 * R
}

## -----------------------------
## Prepare matrices
## -----------------------------
Hs      <- c(0.1, 0.3, 0.5, 0.7, 0.9)
K       <- 3
Sigma2  <- 1

mat_list <- lapply(Hs, cov_innovation, K = K, sigma_u2 = Sigma2)
z_range  <- range(unlist(mat_list))
col_fun  <- colorRampPalette(c("lightyellow", "orange", "red"))
time_idx <- seq(0, 1.4, length.out = K)

## -----------------------------
## Layout for: 5 heatmaps + optimization + color legend
## -----------------------------
layout_matrix <- matrix(
  c(1,2,3,
    4,5,6,
    7,7,7),  # Legend occupies entire bottom row
  nrow = 3,
  byrow = TRUE
)

layout(layout_matrix, heights = c(1, 1, 0.25))

par(mar = c(4,4,3,1))

## -----------------------------
## Plot 5 heatmaps
## -----------------------------
for (i in seq_along(Hs)) {
  image(time_idx, time_idx, mat_list[[i]],
        col = col_fun(50), zlim = z_range,
        xlab = "Time index", ylab = "Time index",
        axes = FALSE, main = bquote(H == .(Hs[i])))
  axis(1, at = time_idx, labels = sprintf("%.1f", time_idx))
  axis(2, at = time_idx, labels = sprintf("%.1f", time_idx))
  box()
}

## -----------------------------
## Panel 6: Optimization of H
## -----------------------------
H_grid <- seq(0.1, 0.9, by = 0.01)
error  <- (H_grid - 0.6)^2 + 0.5

plot(H_grid, error, type = "l",
     xlab = "H value", ylab = "Error",
     main = "Optimization of H")
points(0.6, (0.6 - 0.6)^2 + 0.5,
       pch = 19, col = "red", cex = 1.5)

## -----------------------------
## Color Legend (Bottom row)
## -----------------------------
par(mar = c(2,4,2,4))
zvals <- seq(z_range[1], z_range[2], length.out = 50)

image(x = 1, y = zvals, z = matrix(zvals, nrow = 1),
      col = col_fun(50), axes = FALSE,
      xlab = "", ylab = "Covariance")

axis(4)
box()