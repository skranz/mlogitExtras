#' Get tidy output of estimated coefficients of mlogit model
#'
#' @param mod the estimated model
#' @param scale A numeric value or variable name used to scale coefficients and standard errors. Typical application would be a "price" variable if you estimated a product choice. Then utilities can be interpreted as willingness to pay.
ml_tidy = function(mod, scale=NA) {
  sum = summary(mod)
  co = coef(sum)
  if (!is.na(scale)) {
    if (is.character(scale)) {
      scale = abs(co[scale,1])
    }
    co[,1:2] = co[,1:2] / scale
  }
  as.data.frame(co)
}
