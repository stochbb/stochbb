loadModule("stochbb", TRUE)

is.Var <- function(x) { class(x) == Var; }

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

normal <- function(mu, sigma, name="") {
  if (is.numeric(mu) && is.numeric(sigma)) {
    normalrv(mu, sigma, name);
  } else if (is.Var(mu) && is.numeric(sigma)) {
    compnormalrv(mu, delta(sigma), name);
  } else if (is.numeric(mu) && is.Var(sigma)) {
    compnormalrv(delta(mu), sigma, name);
  } else {
    compnormalrv(mu, sigma, name);
  }
}

gamma <- function(k, theta, name="") {
  if (is.numeric(k) && is.numeric(theta)) {
    gammarv(k, theta, name);
  } else if (is.Var(k) && is.numeric(theta)) {
    compgammarv(k, delta(theta), name);
  } else if (is.numeric(k) && is.Var(theta)) {
    compgammarv(delta(k), theta, name);
  } else {
    compgammarv(k, theta, name);
  }
}

invgamma <- function(alpha, beta, name="") {
  if (is.numeric(alpha) && is.numeric(beta)) {
    invgammarv(alpha, beta, name);
  } else if (is.Var(alpha) && is.numeric(beta)) {
    compinvgammarv(alpha, delta(beta), name);
  } else if (is.numeric(alpha) && is.Var(beta)) {
    compinvgammarv(delta(alpha), beta, name);
  } else {
    compinvgammarv(alpha, beta, name);
  }
}

weibull <- function(k, lambda, name="") {
  if (is.numeric(k) && is.numeric(lambda)) {
    weibullrv(k, lambda, name);
  } else if (is.Var(k) && is.numeric(lambda)) {
    compweibullrv(k, delta(lambda), name);
  } else if (is.numeric(k) && is.Var(lambda)) {
    compweibullrv(delta(k), lambda, name);
  } else {
    compweibullrv(k, lambda, name);
  }
}
