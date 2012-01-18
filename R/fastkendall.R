
## concor calculates all correlations and concordance counts.

concor = function(x, y){
  if(length(x) != length(y))
    stop("x and y must have the same length")
  if(any(is.na(x)) || any(is.na(y)))
    return(NA)
  
  cnt = .Call("concordance_count", as.numeric(x), as.numeric(y))
  result = list(
    concordant=cnt[1], discordant=cnt[2],
    extraX=cnt[3], extraY=cnt[4], spare=cnt[5]
    )

  result$kendall.tau = (result$concordant - result$discordant) / (sqrt(result$concordant + result$discordant + result$extraX) * sqrt(result$concordant + result$discordant + result$extraY))

  ## result$kendall.tau.2 = (result$concordant - result$discordant) / (result$concordant + result$discordant + result$extraX + result$extraY)

  result$goodman.kruskal.gamma = (result$concordant - result$discordant) / (result$concordant + result$discordant)

  ## class(result) = "concor"
  
  return(result)
}
  

concor.benchmark = function(l){
  ## l: vector of vector lenghts

  result = data.frame(l=l, concor=NA, cor.kendall=NA, cor=NA)

  for(i in 1:nrow(result)){
    x = runif(result$l[i])
    y = runif(result$l[i])
    
    result$concor[i] = system.time(concor(x,y))[3]
    result$cor[i] = system.time(cor(x, y))[3]
    result$cor.kendall[i] = system.time(cor(x, y, method="kendall"))[3]
  }

  return(result)
}
