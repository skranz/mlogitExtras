#' Helper function to show changes in choice probabilities / market shares
#'
#' @export
show_P_changes = function(org,new, digits=1) {
  cbind(
    data.frame(var=c("org share","new share","change percentage points", "change percent")),
    round(digits=digits, 100*rbind(
      org,
      new,
      new-org,
      (new-org)/org
    ))
  )
}
