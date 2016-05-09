library(methods)
loadModule("stochbb", TRUE)

is.Var <- function(x) { class(x) == Var; }

uniform <- function(a, b, name="") {
  stochbb::"_uniformrv"(a,b,name);
}

normal <- function(mu, sigma, name="") {
  if (is.numeric(mu) && is.numeric(sigma)) {
    stochbb::"_normalrv"(mu, sigma, name);
  } else if (is.Var(mu) && is.numeric(sigma)) {
    stochbb::"_compnormalrv"(mu, delta(sigma), name);
  } else if (is.numeric(mu) && is.Var(sigma)) {
    stochbb::"_compnormalrv"(delta(mu), sigma, name);
  } else {
    stochbb::"_compnormalrv"(mu, sigma, name);
  }
}

gamma <- function(k, theta, name="") {
  if (is.numeric(k) && is.numeric(theta)) {
    stochbb::"_gammarv"(k, theta, name);
  } else if (is.Var(k) && is.numeric(theta)) {
    stochbb::"_compgammarv"(k, delta(theta), name);
  } else if (is.numeric(k) && is.Var(theta)) {
    stochbb::"_compgammarv"(delta(k), theta, name);
  } else {
    stochbb::"_compgammarv"(k, theta, name);
  }
}

invgamma <- function(alpha, beta, name="") {
  if (is.numeric(alpha) && is.numeric(beta)) {
    stochbb::"_invgammarv"(alpha, beta, name);
  } else if (is.Var(alpha) && is.numeric(beta)) {
    stochbb::"_compinvgammarv"(alpha, delta(beta), name);
  } else if (is.numeric(alpha) && is.Var(beta)) {
    stochbb::"_compinvgammarv"(delta(alpha), beta, name);
  } else {
    stochbb::"_compinvgammarv"(alpha, beta, name);
  }
}

weibull <- function(k, lambda, name="") {
  if (is.numeric(k) && is.numeric(lambda)) {
    stochbb::"_weibullrv"(k, lambda, name);
  } else if (is.Var(k) && is.numeric(lambda)) {
    stochbb::"_compweibullrv"(k, delta(lambda), name);
  } else if (is.numeric(k) && is.Var(lambda)) {
    stochbb::"_compweibullrv"(delta(k), lambda, name);
  } else {
    stochbb::"_compweibullrv"(k, lambda, name);
  }
}

studt <- function(nu, name="") {
  if (is.numeric(k)) {
    stochbb::"_studtrv"(nu, name);
  } else {
    stochbb::"_compstudtrv"(nu, name);
  }
}

chain <- function(X1, X2, ...) {
  args <- c(X1, X2, ...)
  stochbb::"_chain"(args);
}

independent <- function(X1, X2, ...) {
  args <- c(X1, X2, ...)
  stochbb::"_independent"(args);
}

maximum <- function(X1, X2, ...) {
  args <- c(X1, X2, ...)
  stochbb::"_maximum"(args);
}

minimum <- function(X1, X2, ...) {
  args <- c(X1, X2, ...)
  stochbb::"_minimum"(args);
}

'%+%' <- function (X1, X2) {
  if (is.numeric(X1) && is.Var(X2)) {
    affine(X2, 1, X1)
  } else if (is.numeric(X2) && is.Var(X1)) {
    affine(X1, 1, X2)
  } else if (is.Var(X1) && is.Var(X2)) {
    chain(X1,X2)
  } else {
    X1+X2;
  }
}

