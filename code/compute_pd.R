compute_pd <- function(x, a, b, lambda) {
  # This function computes the lambda-PD statistics for a response pattern x at
  # theta = -3.0, -2.99, -2.98, ..., 3.0
  theta <- seq(from = -3, to = 3, by = 0.01)
  p <- 1 / (1 + exp(a * outer(b, theta, "-")))
  e_1 <- apply(a * p, 2, sum)
  e_2 <- apply(a * (1 - p), 2, sum)
  a_sum <- sum(a)
  var_n_1 <- apply(a^2 * p * (1 - p), 2, sum)
  out <- 2 * e_1 * (a_sum - e_1) / (var_n_1 * a_sum * lambda * (lambda + 1)) *
  (n_1 * ((n_1 / e_1)^lambda - 1) + n_2 * ((n_2 / e_2)^lambda - 1))
  return(out)
}
n_item <- 10
a <- rlnorm(n_item, 0, 0.15)
b <- rnorm(n_item, 0, 1)
irt_2pl <- function(theta, a, b) {
  out <- 1 / (1 + exp(a * outer(b, theta, "-")))
  return(out)
}
x <- 1 * (irt_2pl(-1, a, b) > 0.5) # generate a typical response pattern
pd <- compute_pd(x, a, b, lambda = 1)
plot(pd ~ theta, type = "l", ylim = c(-0.5, 1))
abline(h = qchisq(p = 0.025, df = 1))