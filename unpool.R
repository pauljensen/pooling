
require(dplyr)

# convert codes to pooling matrix
codes <- read.csv("codes.csv", colClasses=c("integer", "integer", "character"))
n.wells <- dim(codes)[[1]]
n.pools <- nchar(codes$code[1])
pools <- matrix(FALSE, nrow=n.pools, ncol=n.wells)
for (i in 1:n.wells) {
  pools[,i] <- strsplit(codes$code[i], "")[[1]] == "1"
}
# if pools[i,j] == TRUE, then well j was placed into pool i

pool.names <- c("A1","B1","C1","D1",
                "A2","B2","C2","D2",
                "A3","B3","C3","D3",
                "A4","B4","C4","D4",
                "A5","B5","C5","D5",
                "A6","B6","C6","D6")

load.pool <- function(pool) {
  # pool is a string containing the pool name "A1", "C2", etc.
  # Returns a data frame with the raw read counts.
  filename = paste("reads/Pool_", pool, ".csv", sep="")
  raw <- read.csv(filename)[,c("position","count_1","gene")]
  data.frame(position=raw$position, count=raw[,"count_1"],
             gene=raw$gene, pool=pool)
}

# combine read counts from all pools
reads <- data.frame()
for (pool in pool.names) {
  reads <- rbind(reads, load.pool(pool))
}

# we should ignore insertion sites with fewer than READ.CUTOFF reads
READ.CUTOFF <- 15
reads <- filter(reads, count >= READ.CUTOFF)

# create table mapping positions to genes
positions <- distinct(select(reads, position, gene))  
n.positions <- nrow(positions)

# if trans[i,j] == 1, then transposon insertion site j was present in pool i
trans <- matrix(FALSE, nrow=n.pools, ncol=n.positions)
row.names(trans) <- pool.names
for (i in 1:n.positions) {
  hits <- reads$pool[reads$position == positions$position[i]]
  trans[hits,i] <- TRUE
}
row.names(trans) <- NULL

find.closest <- function(pools, trans, status=F) {
  # Find the pooling pattern that most closely matches the pattern of
  # read abundances.
  #
  # pools is an (n.pools x n.colonies) matrix describing the pooling patterns
  # trans is an (n.pools x n.positions) matrix describing read presence
  # status is a logical controlling status reporting
  #
  # Returns a data frame with n.positions rows.  trans.index is the column
  # index in trans.  pool.index is the closest-matching column index in pools.
  # distance is the Hamming distance between the read presence and the pooling
  # pattern.
  pooling <- t(pools)
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

# this is it; decode all the pools
key <- find.closest(pools, trans, status=T)

# add position and gene lookups
all.key <- cbind(key, positions)
all.key <- cbind(all.key, codes[all.key$pool.index,])
row.names(all.key) <- NULL

# print the number of positions and genes matched at each threshold
for (thresh in 0:4) {
  hits <- filter(all.key, distance <= thresh)
  print(c(nrow(hits %>% select(position) %>% distinct), 
          nrow(hits %>% select(gene) %>% distinct)))
}

# apply the Hamming distance cutoff and output the keys
filtered.keys.3 <- all.key %>% filter(distance <= 3) %>% arrange(plate, well)
write.csv(filtered.keys.3, file="filtered_distance3.csv")
