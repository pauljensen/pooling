

find.closest <- function(pool, trans, status=F) {
  pooling <- t(pool)
  n.trans <- dim(trans)[2]
  key <- data.frame(trans.index=1:n.trans, pool.index=0, distance=0)
  for (i in 1:n.trans) {
    pattern <- as.vector(trans[,i])
    distances <- apply(pooling,1,function(x) sum(!(x == pattern)))
    key$pool.index[i]=which.min(distances)
    key$distance[i]=min(distances)
    
    # status
    if (status) {
      if (i %% 1000 == 0)
        print(i/n.trans * 100)
    }
  }
  key
}
