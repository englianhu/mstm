
# Simulating multivariate aerial spatiotemporal data ----------------------
library(clusterGeneration)
library(spdep)
library(tidyverse)
library(viridis)

# Number of dimensions
L <- 2

# Number of timesteps
Nt <- 10

example(columbus)
coords <- coordinates(columbus)
xx <- poly2nb(columbus)
plot(columbus, border = "grey")
plot(col.gal.nb, coords, add = TRUE)

columbus_map <- fortify(columbus) %>%
  mutate(id = as.numeric(id) + 1)

columbus_df <- as.data.frame(columbus)


# create a multivariate spatiotemporal version of columbus data
mvdf <- expand.grid(POLYID = 1:nrow(columbus_df), l = 1:L, t = 1:Nt)

d <- left_join(mvdf, columbus_df) %>%
  tbl_df %>%
  arrange(t)

# Make some temporally varying covariates
d <- d %>%
  mutate(hoval = HOVAL * sqrt(t) + rnorm(n(), 0, 5), 
         inc = INC * -t^2 + rnorm(n(), 0, 5), 
         plumb = PLUMB * .3 * t + rnorm(n(), 0, 5), 
         choval = c(scale(hoval)), 
         cinc = c(scale(inc)), 
         cplumb = c(scale(plumb)))

pairs(d[, c("choval", "cinc", "cplumb")])


W <- as.matrix(nb2mat(xx, style = "B"))

# Simulate one draw from a CAR distribution
rcar <- function(adjacency_matrix, tau, alpha) {
  # draw sample from proper car model
  # phi ~ N(0, (tau * (D - alpha * W)) ^ -1)
  D <- diag(rowSums(adjacency_matrix))
  Sigma <- solve(tau * (D - alpha * adjacency_matrix))
  L_Sigma <- t(chol(Sigma))
  c(L_Sigma %*% rnorm(nrow(adjacency_matrix)))
}

x <- rcar(W, 1, .99)

columbus_map %>%
  mutate(zval = x[id]) %>%
  ggplot(aes(x = long, y = lat, group = group, fill = zval)) + 
  geom_polygon() + 
  scale_fill_viridis()

# Now put together design matrices for each timestep
X <- array(dim = c(nrow(d) / Nt, 2, Nt))
for (i in 1:Nt) {
  X[, , i] <- d %>%
    filter(t == i) %>%
    model.matrix(~ 1 + cplumb, data = .)
}

# Make block diagonal adjacency matrix
A_t <- matrix(0, nrow = nrow(X[, , 1]), ncol = nrow(X[, , 1]))
A_t[1:nrow(columbus_df), 1:nrow(columbus_df)] <- W
A_t[(nrow(columbus_df) + 1):nrow(A_t), (nrow(columbus_df) + 1):nrow(A_t)] <- W
image(A_t)

# compute multivariate basis functions
r <- 5
S_X <- array(dim = c(nrow(X[,,1]), r, Nt))
for (i in 1:Nt) {
  G <- (diag(nrow(X[, , 1])) - 
          X[, , i] %*% solve(t(X[, , i]) %*% X[, , i]) %*% t(X[, , i])) %*%
    A_t %*%
    (diag(nrow(X[, , 1])) - 
       X[, , i] %*% solve(t(X[, , i]) %*% X[, , i]) %*% t(X[, , i]))
  eG <- eigen(G)
  basis_vectors <- eG$vectors
  S_X[, , i] <- basis_vectors[, 1:r]
}

# visualize multivariate basis functions
G_unrolled <- S_X[, , 1]
for (i in 2:Nt) G_unrolled = rbind(G_unrolled, S_X[, , i])

tbl_df(G_unrolled) %>%
  bind_cols(d) %>%
  mutate(id = POLYID) %>%
  right_join(columbus_map) %>%
  ggplot(aes(x = long, y = lat, group = id, fill = V4)) + 
  geom_polygon() + 
  facet_grid(l ~ t) + 
  scale_fill_viridis() + 
  coord_equal()

# confirm that the basis functions are independent of the covariates
tbl_df(G_unrolled) %>%
  bind_cols(d) %>%
  select(starts_with("V"), choval, cinc, cplumb) %>%
  pairs

# construct the propagator matrix
Bt <- array(dim = c(r, dim(X)[2] + r, Nt))
Mt <- array(dim = c(r, r, Nt))
for (i in 1:Nt) {
  Bt[, , i] <- cbind(t(S_X[, , i]) %*% X[, , i], diag(r))
  G_B <- (diag(r) - 
            Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i])) %*%
    diag(r) %*% 
    (diag(r) - 
       Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i]))
  eGB <- eigen(G_B)
  Mt[, , i] <- eGB$vectors[, 1:r]
}

# construct the basis function coefficients
eta <- matrix(nrow = Nt, ncol = r)
R_eta <- rcorrmatrix(r, alphad = 5)
LR_eta <- t(chol(R_eta))
sigma_eta <- .1
eta[1, ] <- rnorm(r)
for (i in 2:Nt) {
  eta[i, ] <- Mt[, , i] %*% eta[i - 1, ] + sigma_eta * c(LR_eta %*% rnorm(r))
}

# coef for fixed effects
beta <- rnorm(dim(X)[2])

# process model
Y <- matrix(nrow = dim(X)[1], ncol = Nt)
for (i in 1:Nt) {
  Y[, i] <- X[, , i] %*% beta + S_X[, , i] %*% eta[i, ]
}

# observation model
Z <- Y + matrix(rnorm(dim(X)[1] * Nt, sd = .2), ncol = Nt)
