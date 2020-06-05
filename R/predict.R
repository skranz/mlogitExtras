
#' Alternative predict function for mlogit objects
#'
#' Currently only works for models without alternative specific constants or
#' alternative specific interaction effects.
#'
#' @param mod An estimated mlogit model
#' @param newdata A data set for the prediction in long format. The data set should have a column alt the indexed the alternative. The data set must be ordered first by choice situation and each choice situation must have the same number of alternatives in the same order. You may simply add a column chid specifying choice situation then arrange by chid, alt.
#' @param num.draws Number of simulated consumers to compute market shares (relevant for mixed logit models only)
#' @param use.halton Should Halton sequences instead of pseudo-random numbers be used to simulate consumers (default TRUE)
#' @return A matrix of predicted choice probabilities with one row per choice situation and one column for each alternative. Each row sums up to 1.
ml_predict = function(mod, newdata,num.draws = 1000, use.halton = TRUE) {
  restore.point("mlogit_predict")

  coef = mod$coefficients
  coef = coef[!startsWith(names(coef), "sd.")]
  rvars = names(mod$rpar)

  if (is.null(newdata[["alt"]])) {
    stop("Your data set must have a column alt that specify the choice alternative. It must be ordered by choice situations and each choice situations must have the same alternatives in the same order.")
  }

  alts = unique(newdata$alt)
  J = length(alts)

  # This only works for mlogit models without alternative specific
  # constant or any other alternative specific interaction effects
  X = model.matrix(mod$formula, newdata)
  if (colnames(X)[1]== "(Intercept)")
    X = X[,-1]

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
  fvars = fvars[!startsWith(fvars,"chol.") & !startsWith(fvars,"sd.")]

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
