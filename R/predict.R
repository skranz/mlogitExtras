
#' Alternative predict function for mlogit objects
ml_predict = function(mod, newdata,J=NA, num.draws = 1000, alts=NA, use.halton = TRUE) {
  restore.point("mlogit_predict")

  coef = mod$coefficients
  coef = coef[!startsWith(names(coef), "sd.")]
  rvars = names(mod$rpar)

  if (is.na(J)) {
    if (is.null(newdata$chid)) {
      stop("Either proivde the number of alternatives J or have a column chid.")
    }
    J = sum(newdata$chid == newdata$chid[[1]])
  }

  if (all(is.na(alts))) {
    if (!is.null(newdata$alt)) {
      alts = newdata$alt[1:J]
    } else {
      alts = NULL
    }
  }

  #X = model.matrix(mod$formula, newdata)
  #if (colnames(X)[1]== "(Intercept)")
  #  X = X[,-1]

  X = as.matrix(newdata[, names(coef)])
  # No random coefficients
  if (length(rvars)==0) {
    V = matrix(X %*% coef, ncol=J, byrow=TRUE)
    expV = exp(V)
    rowSumsExpV = rowSums(expV)
    P = expV / rowSumsExpV
    colnames(P) = alts
    return(P)

  }
  # We have random coefficients

  num.sit = NROW(newdata) / J

  # Non-random variables
  fvars = setdiff(names(coef),rvars)

  # Draw individual betas
  beta.mat = ml_draw_betas(mod, num.draws=num.draws, use.halton = use.halton)

  # Resulting probability matrix
  P = matrix(NA, nrow=num.sit, ncol=J)
  colnames(P) = alts
  sit = 1
  for (sit in 1:num.sit) {
    rows = (sit-1)*J+(1:J)
    if (length(fvars)>0) {
      V = as.vector(X[rows,fvars,drop=FALSE] %*% coef[fvars])
    } else {
      V = numeric(J)
    }
    V.mat = matrix(V, nrow=num.draws, ncol=J, byrow = TRUE)
    k = 1
    for (k in seq_along(rvars)) {
      x = X[rows,rvars[k]]
      x.temp = matrix(x, nrow=num.draws, ncol=J,byrow=TRUE)
      V.mat = V.mat + x.temp * beta.mat[,k]
    }
    expV = exp(V.mat)
    rowSumsExpV = rowSums(expV)
    P.mat = expV / rowSumsExpV

    P[sit,] = colMeans(P.mat)
  }
  P
}
