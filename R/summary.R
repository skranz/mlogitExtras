#' Get tidy output of estimated coefficients of mlogit model as data frame
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
  vars = rownames(co)
  is.chol = startsWith(vars, "chol.")
  if (any(is.chol)) {
    co = co[!is.chol,]
    sr = summary(vcov(mod, what="rpar", type="cor"))
    att = attributes(sr)
    sr.mat = matrix(as.numeric(sr),nrow=att$dim[1], ncol=att$dim[2],dimnames = att$dimnames )
    if (!is.na(scale)) {
      sd.rows = startsWith(rownames(sr.mat),"sd.")
      sr.mat[sd.rows,1:2] = sr.mat[sd.rows,1:2] / scale
    }
    co = rbind(co, sr.mat)
  }

  vars = rownames(co)

  rpar = mod[["rpar"]]
  if (length(rpar)>0) {
    dist = sapply(mod[["rpar"]], function(rp) rp$dist)
    inds = startsWith(vars, "sd.")
    co[inds,c(1,3)] = abs(co[inds, c(1,3)])

    span.vars = names(dist)[dist %in% c("u","t")]
    inds = match(paste0("sd.",span.vars), vars)
    vars[inds] = paste0("span.", span.vars)
  }

  rownames(co) = NULL
  colnames(co) = c("estimate","std.error","z.value","p.value")
  p.value = co[,4]
  #sig = ifelse(p.value <= 0.001, "`***`",
#    ifelse(p.value <= 0.01, "`**`",
#    ifelse(p.value <= 0.05, "`*`",
#    ifelse(p.value <= 0.1, ".",
#    ""
#  ))))
  #df = cbind(data.frame(variable = vars),as.data.frame(co), data.frame(signif=sig))
  df = cbind(data.frame(variable = vars),as.data.frame(co))

  df
}
