
library('Matrix')
source("bitmask.R")

set.seed(1986)

pick.wells <- function(n,freq) {
  # n is the number of wells to pick
  # freq is the frequency of each transposon
  sample.int(length(freq),size=n,replace=T,prob=freq)
}

make.pools <- function(wells,n.trans,n.pools) {
  n.wells <- length(wells)
  pool.size <- ceiling(n.wells/2)
  pools <- matrix(F,nrow=n.pools,ncol=n.wells)
  contents <- matrix(F,nrow=n.pools,ncol=n.trans)
  for (i in 1:n.pools) {
    wells.subset <- sample.int(n.wells,size=pool.size,replace=F)
    pools[i,wells.subset] <- T
    contents[i,wells[wells.subset]] <- T
  }
  list(pools=pools, contents=contents)
}

make.crt.pools <- function(wells,n.pools,n.trans) {
  n.windows <- length(n.pools)
  n.wells <- length(wells)
  windows <- vector('list',n.windows)
  pool.size <- ceiling(n.wells/min(n.pools))
  for (i in 1:n.windows) {
    n.pool <- n.pools[i]
    aug.wells <- integer(n.pool * pool.size)
    aug.wells[1:length(wells)] <- wells
    windows[[i]] <- matrix(aug.wells,nrow=n.pool,ncol=pool.size)
  }
  pool.keys <- Reduce(rbind,windows)
  pools <- matrix(F,nrow=sum(n.pools),ncol=n.wells)
  contents <- matrix(F,nrow=sum(n.pools),ncol=n.trans)
  for (i in 1:dim(pool.keys)[[1]]) {
    pools[i,pool.keys[i,]] <- T
    contents[i,wells[pools[i,]]] <- T
  }
  return(list(pools=pools,contents=contents))
}

make.contents <- function(pools,wells,n.trans) {
  n.pools <- dim(pools)[[1]]
  contents <- matrix(F,nrow=n.pools,ncol=n.trans)
  for (i in 1:n.pools)
    contents[i,wells[pools[i,]]] <- T
  return(contents)
}

unique.map.to <- function(mask1,mask2) {
  for.idx <- match.masks(mask1,mask2)
  rev.idx <- match.masks(mask1,mask2,reverse=T)
  for.idx[for.idx != rev.idx] <- NA
  return(for.idx)
}

resolve.pools <- function(pools,contents) {
  pool.mask <- logical.to.mask.matrix(pools)
  contents.mask <- logical.to.mask.matrix(contents)
  unique.map.to(contents.mask,pool.mask)
}

#freq <- map$reads / sum(map$reads)
#n.trans <- length(freq)
#wells <- pick.wells(11000,freq)

#pool.results <- make.pools(wells,n.trans,n.pools=24)
#pool.results <- make.crt.pools(wells,n.pools=c(11,13),n.trans)
#pools <- pool.results[[1]]
#contents <- pool.results[[2]]
#pools <- rbind(Ul[,1:length(wells)],Ul[,2:(length(wells)+1)])
#pools <- Ul[,1:length(wells)]
#contents <- make.contents(pools,wells,n.trans)
#well.key <- resolve.pools(pools,contents)

#print("Unique genes in wells")
#print(length(unique(wells)))
#print("Unambiguous wells identified")
#print(sum(!is.na(well.key)))
