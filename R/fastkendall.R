
## fastkendall calculates all correlations and concordance counts.

fastkendall = function(x, y){
  if(length(x) != length(y))
    stop("x and y must have the same length")
  if(any(is.na(x)) || any(is.na(y)))
    return(NA)
  if(!is.numeric(x) || !is.numeric(y))
    stop("x and y must be numeric")
  
  cnt = .Call("concordance_count", as.numeric(x), as.numeric(y))
  result = list(
    concordant=cnt[1], discordant=cnt[2],
    extra.x=cnt[3], extra.y=cnt[4], spare=cnt[5]
    )

  result$kendall.tau = (result$concordant - result$discordant) / (sqrt(result$concordant + result$discordant + result$extra.x) * sqrt(result$concordant + result$discordant + result$extra.y))

  ## result$kendall.tau.2 = (result$concordant - result$discordant) / (result$concordant + result$discordant + result$extra.x + result$extra.y)

  result$goodman.kruskal.gamma = (result$concordant - result$discordant) / (result$concordant + result$discordant)

  return(result)
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
