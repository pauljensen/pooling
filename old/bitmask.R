
logical.to.mask.matrix <- function(m, max.bits=28) {
  r <- dim(m)[[1]]
  c <- dim(m)[[2]]
  if (r <= max.bits) {
    scales <- 2^(0:(r-1))
    return(as.integer(apply(m,2,function(x) sum(x * scales))))
  }
  
  # make the matrix an even multiple of max.bits
  n.div <- ceiling(r/max.bits)
  new.r <- n.div * max.bits
  m <- rbind(m, matrix(F,nrow=new.r-r,ncol=c))
  
  masks <- matrix(0L,nrow=n.div,ncol=c)
  scales <- 2^(0:(max.bits-1))
  starts <- seq(1,new.r,max.bits)
  for (i in 1:n.div) {
    masks[i,] <- as.integer(apply(m[starts[i]:(starts[i]+max.bits-1),],2,
                                  function(x) sum(x * scales)))
  }
  return(masks)
}

match.masks <- function(mask1,mask2,reverse=F) {
  if (reverse) {
    if (is.vector(mask2)) {
      nc <- length(mask2)
      return(nc + 1 - match.masks(mask1,mask2[seq(nc,1,-1)]))
    } else {
      # matrix
      nc <- dim(mask2)[[2]]
      return(nc + 1 - match.masks(mask1,mask2[,seq(nc,1,-1)]))
    }
  }
  
  if (is.vector(mask1) && is.vector(mask2))
    return(match(mask1,mask2))
  else if (is.vector(mask1))
    mask1 <- matrix(mask1,nrow=1)
  else if (is.vector(mask2))
    mask2 <- matrix(mask2,nrow=1)
  
  key1 <- as.vector(mask1[1,])
  key2 <- as.vector(mask2[1,])
  ord <- order(key2)
  key2s <- key2[ord]
  mask2s <- mask2[,ord]
  first <- match(key1,key2s,nomatch=0)
  for (i in 1:length(first)) {
    if (first[i] == 0) {
      next
    }
    if (all(mask1[,i] == mask2s[,first[i]])) {
      next
    } else {
      # look for other matches with common first word
      k <- first[i] + 1
      found.match <- F
      while (k <= length(key2s) && key2s[k] == key2s[first[i]]) {
        if (all(mask1[,i] == mask2s[,k])) {
          found.match <- T
          break
        }
        k <- k + 1
      }
      if (found.match) {
        first[i] <- k
      } else {
        first[i] <- 0
      }
    }
  }
  
  # un-sort the indices
  first[first>0] <- ord[first[first>0]]
  return(first)
}

set.seed(1988)
cont <- matrix(runif(30)>0.5, nrow=6, ncol=5)

m1 <- logical.to.mask.matrix(cont)
m2 <- logical.to.mask.matrix(cont,max.bits=2)

print(m2)
#print(m2[,-1])
print(match.masks(m2,m2[,-2]))

