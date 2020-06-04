#' Draw random coefficients of individuals for an estimated mixed logit model
#'
#' Note that the function is still restricted and does not yet work for all sorts
#' of random parameter distributions that mlogit supports.
#'
#' @param mod the estimated models
#' @param num.draws the number of individuals for which random coefficients shall be drawn
#' @param use.halton shall Halton sequences be used instead of psedudo-random numbers?
#' @return A matrix with num.draws rows and one column for each explanatory variable that has random coefficients
ml_draw_betas = function(mod, num.draws=100, use.halton = FALSE, scale=NA) {
  restore.point("ml_draw_betas")

  rpar=mod$rpar
  rvars = names(rpar)
  K = length(rvars)

  chol = ml_chol(mod)
  cvars = colnames(chol)
  ivars = setdiff(rvars, cvars)

  beta.mat = matrix(NA, nrow=num.draws, ncol=K)
  colnames(beta.mat) = rvars

  if (use.halton) {
    beta.unif = randtoolbox::halton(num.draws, dim=K)
  } else {
    beta.unif = matrix(runif(num.draws*K), nrow=num.draws, ncol=K)
  }
  colnames(beta.unif) = rvars

  if (length(cvars)>0) {
    dist = sapply(rpar[cvars], function(rp) rp$dist)

    if (any(!dist %in% c("n","cn"))) {
      stop("The mlogitExtras package can currently only deal with correlated betas that are normally distributed or zero censored normally distributed.")
    }

    # First convert uniform betas into standard normal betas
    beta.norm = qnorm(beta.unif[,cvars])

    # Use cholesky matrix to create correlated normal variables
    beta.norm = t(chol %*% t(beta.norm))

    # Add means of betas
    means = sapply(rpar[cvars], function(rp) rp$mean)
    beta.mat[, cvars] = beta.norm + matrix(means,num.draws, length(cvars), byrow = TRUE)

    if (any(dist=="cn")) {
      cols = cvars[dist == "cn"]
      inds = beta.mat[,cols] < 0
      beta.mat[inds] = 0
    }

    # mat = t(chol %*% t(beta.norm))
    # var(mat)
    # cor(mat)
    # Kc = length(cvars)
    #
    # summary(mod)
    # rpar(mod)
    # vcov(mod, what="rpar")

  }

  # Transform remaining iid betas
  for (ivar in ivars) {
    beta.mat[, ivar] = transform.uniform.to.rpar(beta.unif[,ivar],mod$rpar[[ivar]])
  }

  # Possibly scale the betas
  if (!is.na(scale)) {
    if (is.character(scale)) {
      scale = abs(coef(mod)[scale])
    }
    beta.mat = beta.mat / scale
  }
  beta.mat
}

draw.iid.rpar.beta = function(rpar, num.draws=1000) {
  if (rpar$dist == "n") {
    rnorm(num.draws,rpar$mean, abs(rpar$sigma))
  } else if (rpar$dist == "u") {
    half.span = rpar$sigma / 2
    runif(num.draws, rpar$mean-half.span, rpar$mean+half.span)
  } else {
    stop("Prediction is implemented so far only for iid normal and iid uniform distributed random coefficients.")
  }
}

# Transform uniformely distributed variables
transform.uniform.to.rpar = function(beta.unif, rpar) {
  if (rpar$dist == "n") {
    qnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
  } else if (rpar$dist == "u") {
    span = rpar$sigma
    beta.unif*span+(rpar$mean-span/2)
  } else if (rpar$dist == "t") {
    min = rpar$mean-0.5*abs(rpar$sigma)
    max = rpar$mean+0.5*abs(rpar$sigma)
    qtri(beta.unif, min,max)
  } else if (rpar$dist == "ln") {
    qlnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
  } else if (rpar$dist == "cn") {
    beta = qnorm(beta.unif, mean = rpar$mean, sd = abs(rpar$sigma))
    beta[beta<0] = 0
    beta
  } else {
    stop("Prediction is implemented so far only for iid normal and iid uniform distributed random coefficients.")
  }


}

#' Extract Cholesky decomposition matrix of correlated random variables of an mlogit model
ml_chol = function(mod) {
  co = coef(mod)
  co = co[startsWith(names(co),"chol.")]
  if (length(co)==0) return(NULL)
  na = substring(names(co),6)
  var.mat = matrix(unlist(strsplit(na, ":", fixed=TRUE)), ncol=2, byrow=TRUE)
  var.mat
  vars = unique(var.mat[,1])
  num.var = length(vars)
  chol = matrix(0,num.var, num.var)
  rownames(chol) = colnames(chol) = vars
  chol[var.mat] = co
  t(chol)
}

