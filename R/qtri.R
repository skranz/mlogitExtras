example.qtri = function() {
  p = seq(0,1,by=0.001)
  q = qtri(p,max=2)
  d = dtri(q,max=2)
  cumsum(d[-1]*diff(q))
  plot(x=q,y=d)

}

#' Density function of triangular distribution
dtri <- function (x, min = 0, max = 1){
  n = max(length(x),length(min), length(max))
  if (n>length(x)) p = rep(x,length.out=n)
  if (n>length(min)) min = rep(min,length.out=n)
  if (n>length(max)) max = rep(max,length.out=n)

  names.x <- names(x)
  d = numeric(n)
  if (any(x < min) || any(x > max))
      stop("All values of 'x' must be between min and max.")
  if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All values of 'min' and 'max' must be finite.")

  mode = (max+min)/2
  range = (max-min)
  d <- 2 * ifelse(x <= mode,
    (x - min)/(range * (mode -min)),
    (max - x)/(range * (max - mode))
  )
  if (!is.null(names.x))
    names(d) <- rep(names.x, length = length(d))
  else names(d) <- NULL

  d
}

#' Quantile function of triangular distribution
qtri <- function (p, min = 0, max = 1)  {
  n = max(length(p),length(min), length(max))
  if (n>length(p)) p = rep(p,length.out=n)
  if (n>length(min)) min = rep(min,length.out=n)
  if (n>length(max)) max = rep(max,length.out=n)

  names.p <- names(p)
  q = numeric(n)
  if (any(p < 0) || any(p > 1))
      stop("All values of 'p' must be between 0 and 1.")
  if (any(is.infinite(min)) || any(is.infinite(max)))
      stop("All values of 'min' and 'max' must be finite.")

  q[p == 0] <- min[p == 0]
  q[p == 1] <- max[p == 1]

  range = max-min
  mode = (max+min) / 2

  rows = p<=1/2 & p > 0
  q[rows] = min[rows] + sqrt(range[rows] * (mode[rows]-min[rows]) * p[rows])
  rows = p>1/2 & p < 1
  q[rows] = max[rows] - sqrt(range[rows] * (max[rows]-mode[rows]) * (1-p[rows]))
  q
  if (!is.null(names.p))
    names(q) <- rep(names.p, length = length(q))
  else names(q) <- NULL
  q
}
