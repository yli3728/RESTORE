# The same model is used for every dataset. Only lambda and numerical
# hyperparameters differ between configurations.
#
#   kappa_ij(lambda) = kappa0 * ((1 - lambda) + lambda * w_ij)
#
# lambda = 1.0: fully adaptive Z penalty (lung)
# lambda = 0.5: half fixed / half adaptive Z penalty (brain)

kappa_lambda <- function(kappa0, lambda, adaptive_weight) {
  stopifnot(length(kappa0) == 1L, kappa0 > 0,
            length(lambda) == 1L, lambda >= 0, lambda <= 1,
            all(is.finite(adaptive_weight)), all(adaptive_weight >= 0))
  kappa0 * ((1 - lambda) + lambda * adaptive_weight)
}

sigmoid <- function(x) 1 / (1 + exp(-x))

