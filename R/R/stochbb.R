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

chain <- function(a, b, ...) {
  args <- c(a, b, ...)
  stochbb::"_chain"(args);
}

independent <- function(a, b, ...) {
  args <- c(a, b, ...)
  stochbb::"_independent"(args);
}

maximum <- function(a, b, ...) {
  args <- c(a, b, ...)
  stochbb::"_maximum"(args);
}

minimum <- function(a, b, ...) {
  args <- c(a, b, ...)
  stochbb::"_minimum"(args);
}

'%+%' <- function (a, b) {
  if (is.numeric(a) && (class(b)==Var)) {
    affine(b, 1, a)
  } else if (is.numeric(b) && (class(a)==Var)) {
    affine(a, 1, b)
  } else if ((class(b)==Var) && (class(b)==Var)) {
    chain(a,b)
  } else {
    a+b;
  }
}

