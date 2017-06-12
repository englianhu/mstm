
# Simulating multivariate aerial spatiotemporal data ----------------------
# based on Bradley, Holan, and Wikle 2015
# Multivariate spatio-temporal models for high-dimensional areal data with 
#   application to Longitudinal Employer-Household Dynamics
# https://arxiv.org/abs/1503.00982

library(reshape2)
library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

library(clusterGeneration)
library(spdep)
library(tidyverse)
library(viridis)

# Number of dimensions
L <- 2

# Number of timesteps
Nt <- 20

# Load "Columbus" shapefile and compute neighbors
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
  mutate(hoval = HOVAL, 
         inc = INC, 
         plumb = PLUMB, 
         choval = c(scale(scale(hoval) * scale(t)^2 + plogis(t))), 
         cinc = c(scale(inc * -t + rnorm(n()))), 
         cplumb = c(scale(scale(plumb) * plogis(t) + rnorm(n(), sd = .1))))

d %>%
  select(POLYID, starts_with("c", ignore.case = FALSE), t) %>%
  gather(var, value, -POLYID, -t) %>%
  ggplot(aes(t, value, group = POLYID)) + 
  facet_wrap(~var) + 
  geom_line()

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
  scale_fill_viridis("CAR random effect") +
  coord_equal()

# Now put together design matrices for each timestep
X <- array(dim = c(Nt, nrow(d) / Nt, 2))
for (i in 1:Nt) {
  X[i, , ] <- d %>%
    filter(t == i) %>%
    model.matrix(~ 1 + cinc, data = .)
}

# Make block diagonal adjacency matrix
A_t <- bdiag(replicate(L, list(W))) %>%
  as.matrix
image(A_t)

# compute multivariate basis functions
r <- 5
S_X <- array(dim = c(Nt, nrow(d) / Nt, r))
for (i in 1:Nt) {
  G <- (diag(nrow(X[1, , ])) - 
          X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ])) %*%
    A_t %*%
    (diag(nrow(X[1, , ])) - 
       X[i, , ] %*% solve(t(X[i, , ]) %*% X[i, , ]) %*% t(X[i, , ]))
  eG <- eigen(G)
  basis_vectors <- eG$vectors
  S_X[i, , ] <- basis_vectors[, 1:r]
}

# visualize multivariate basis functions
G_unrolled <- S_X[1, , ]
for (i in 2:Nt) G_unrolled = rbind(G_unrolled, S_X[i, , ])

# V1 is the first basis vector, V2 the second, etc.
tbl_df(G_unrolled) %>%
  bind_cols(d) %>%
  mutate(id = POLYID) %>%
  right_join(columbus_map) %>%
  ggplot(aes(x = long, y = lat, group = id, fill = V1)) +
  geom_polygon() +
  facet_grid(l ~ t) +
  scale_fill_viridis() +
  coord_equal()

# confirm that the basis functions are independent of the covariates
tbl_df(G_unrolled) %>%
  bind_cols(d) %>%
  select(starts_with("V"), cinc) %>%
  cor %>%
  image

# construct the propagator matrix
Bt <- array(dim = c(r, dim(X)[3] + r, Nt))
Mt <- array(dim = c(Nt, r, r))
for (i in 1:Nt) {
  Bt[, , i] <- cbind(t(S_X[i, , ]) %*% X[i, , ], diag(r))
  G_B <- (diag(r) - 
            Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i])) %*%
    diag(r) %*% 
    (diag(r) - 
       Bt[, , i] %*% MASS::ginv(t(Bt[, , i]) %*% Bt[, , i]) %*% t(Bt[, , i]))
  eGB <- eigen(G_B)
  Mt[i, , ] <- eGB$vectors[, 1:r]
}

# simulate basis function coefficients
eta <- matrix(nrow = Nt, ncol = r)
R_eta <- rcorrmatrix(r, alphad = 1)
LR_eta <- t(chol(R_eta))
sigma_eta <- 2
sigma_0 <- 4
eta[1, ] <- sigma_0 * c(LR_eta %*% rnorm(r))
for (i in 2:Nt) {
  eta[i, ] <- Mt[i, , ] %*% eta[i - 1, ] + sigma_eta * c(LR_eta %*% rnorm(r))
}

# simulate coef for fixed effects
(beta <- rnorm(dim(X)[3]))

# process model
Y <- matrix(nrow = dim(X)[2], ncol = Nt)
for (i in 1:Nt) {
  Y[, i] <- X[i, , ] %*% beta + S_X[i, , ] %*% eta[i, ]
}

# visualize latent process
Y %>%
  tbl_df %>%
  gather(year, y) %>%
  bind_cols(d) %>%
  mutate(id = POLYID) %>%
  right_join(columbus_map) %>%
  ggplot(aes(x = long, y = lat, group = id, fill = y)) +
  geom_polygon() +
  facet_grid(l ~ t) +
  scale_fill_viridis() +
  coord_equal() +
  theme_minimal()



# observation model (process + noise)
Z <- Y + matrix(rnorm(dim(X)[2] * Nt, sd = .2), ncol = Nt)

Z %>%
  tbl_df %>%
  gather(year, z) %>%
  bind_cols(d) %>%
  mutate(id = POLYID) %>%
  right_join(columbus_map) %>%
  ggplot(aes(x = long, y = lat, group = id, fill = z)) +
  geom_polygon() +
  facet_grid(l ~ t) +
  scale_fill_viridis() +
  coord_equal()



# Fit model in Stan, evaluate parameter recovery --------------------------

stan_d <- list(n = nrow(columbus_df), 
               T = Nt, 
               L = L, 
               p = dim(X)[3], 
               r = r, 
               S = S_X, 
               M = Mt, 
               Z = Z)

m_init <- stan_model("stan/mstm.stan")

m_fit <- sampling(m_init, data = stan_d)

# evaluate convergence and check recovery of fixed effects
traceplot(m_fit)
beta


post <- rstan::extract(m_fit)


# Evaluate recovery of spatiotemporal basis coefficients
eta_df <- eta %>%
  reshape2::melt(varnames = c("year", "basis_dim"), 
                 value.name = "true_eta")


post$eta %>%
  reshape2::melt(varnames = c("iter", "basis_dim", "year")) %>%
  tbl_df %>%
  group_by(basis_dim, year) %>%
  summarize(med = median(value), 
            lo = quantile(value, .025), 
            hi = quantile(value, .975)) %>%
  ungroup() %>%
  full_join(eta_df) %>% 
  ggplot(aes(true_eta, med)) + 
  geom_point() + 
  geom_segment(aes(xend = true_eta, y = lo, yend = hi)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  xlab("True basis coefficient") + 
  ylab("Estimated basis coefficient")


# Evaluate recovery of basis coefficient hypoerparameters
plot(m_fit, pars = "sigma_eta") + 
  geom_vline(xintercept = sigma_eta, linetype = "dashed")
plot(m_fit, pars = "sigma_0") + 
  geom_vline(xintercept = sigma_0, linetype = "dashed") + 
  xlim(0, 5)

R_eta_df <- R_eta %>%
  reshape2::melt(varnames = c("row", "col"), value.name = "trueR")

post$R_eta %>%
  reshape2::melt(varnames = c("iter", "row", "col")) %>%
  tbl_df %>%
  group_by(row, col) %>%
  summarize(med = median(value), 
            lo = quantile(value, .025), 
            hi = quantile(value, .975)) %>%
  ungroup() %>%
  full_join(R_eta_df) %>% 
  ggplot(aes(trueR, med)) + 
  geom_point() + 
  geom_segment(aes(xend = trueR, y = lo, yend = hi)) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
  xlab("True basis correlation") + 
  ylab("Estimated basis correlation")
