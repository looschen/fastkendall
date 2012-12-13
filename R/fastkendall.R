fastkendall = function(x, y=NULL){
  ## fastkendall calculates both correlations and concordance counts.
  ## x, y may be matrices or vectors as in R's cor.
  ## In case of a matrix, a list of matrices is returned.

  ## no NA-handling
  if(any(is.na(x)) || any(is.na(y)))
    return(NA)

  ## Error handling copied from R's cor
  if (is.data.frame(y))
    y <- as.matrix(y)
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y))
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x)))
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y)))
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }

  fkendall = function(x, y){
    ## basic correlation of two vectors
    cnt = .Call("concordance_count", as.numeric(x), as.numeric(y))
    result = list(
      concordant=cnt[1], discordant=cnt[2],
      extra.x=cnt[3], extra.y=cnt[4], spare=cnt[5]
      )
    result$kendall.tau = (result$concordant - result$discordant) / (sqrt(result$concordant + result$discordant + result$extra.x) * sqrt(result$concordant + result$discordant + result$extra.y))
    ## result$kendall.tau.2 = (result$concordant - result$discordant) / (result$concordant + result$discordant + result$extra.x + result$extra.y)
    result$goodman.kruskal.gamma = (result$concordant - result$discordant) / (result$concordant + result$discordant)
    result
  }

  prepare.matrix.result = function(n, m, dimnames){
    ## return a list of n*m matrices
    M = matrix(nrow=n, ncol=m, dimnames=dimnames)
    list(concordant = M, discordant = M, extra.x = M, extra.y = M, spare = M, kendall.tau = M, goodman.kruskal.gamma = M)
  }


  if(is.null(y)){
    ## correlate cols of x
    ncx = ncol(x)
    mresult = prepare.matrix.result(ncx, ncx, list(colnames(x), colnames(x)))
    mnames = names(mresult)
    for(i in seq_len(ncx)){
      for(j in seq_len(i)){
        tmp = fkendall(x[,i], x[,j])
        for(n in mnames){
          mresult[[n]][i, j] = tmp[[n]]
        }
      }
    }

    ## fill in missing values
    utri = upper.tri(mresult$concordant)
    for(i in mnames)
      mresult[[i]][utri] = t(mresult[[i]])[utri]
    ## extra.x and extra.y are not symmetric, but transposes of each other
    mresult$extra.x[utri] = t(mresult$extra.y)[utri]
    mresult$extra.y[utri] = t(mresult$extra.x)[utri]

    return(mresult)
  }else if(is.matrix(x) || is.matrix(y)){
    ## correlate cols of two matrices
    if(! is.matrix(x))
      x = as.matrix(x)
    if(! is.matrix(y))
      y = as.matrix(y)
    if(nrow(x) != nrow(y))
      stop("'x' and 'y' must have equal number of rows")

    ncx = ncol(x)
    ncy = ncol(y)
    mresult = prepare.matrix.result(ncx, ncy, list(colnames(x), colnames(y)))
    mnames = names(mresult)
    for(i in seq_len(ncx)){
      for(j in seq_len(ncy)){
        tmp = fkendall(x[,i], y[,j])
        for(n in mnames){
          mresult[[n]][i, j] = tmp[[n]]
        }
      }
    }
    return(mresult)
  }else{
    if(length(x) != length(y))
      stop("'x' and 'y' must have equal length")
    return(fkendall(x, y))
  }

}


## fastkendall.benchmark = function(l){
##   ## l: vector of vector lengths

##   result = data.frame(l=l, fastkendall=NA, cor.kendall=NA, cor=NA)

##   for(i in 1:nrow(result)){
##     x = runif(result$l[i])
##     y = runif(result$l[i])

##     result$fastkendall[i] = system.time(fastkendall(x,y))[3]
##     result$cor[i] = system.time(cor(x, y))[3]
##     result$cor.kendall[i] = system.time(cor(x, y, method="kendall"))[3]
##   }

##   return(result)
## }
