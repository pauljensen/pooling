
require(dplyr)

# convert codes to pooling matrix
codes <- read.csv("codes.csv", colClasses=c("integer", "integer", "character"))
n.wells <- dim(codes)[[1]]
n.pools <- nchar(codes$code[1])
pools <- matrix(F, nrow=n.pools, ncol=n.wells)

for (i in 1:n.wells) {
  pools[,i] <- strsplit(codes$code[i], "")[[1]] == "1"
}
# if pools[i,j] == TRUE, then well j was placed into pool i

pool.names <- c("A1","A2","A3","A4","A5","A6",
                "B1","B2","B3","B4","B5","B6",
                "C1","C2","C3","C4","C5","C6",
                "D1","D2","D3","D4","D5","D6")

load.pool <- function(pool) {
  filename = paste("reads/Pool_", pool, ".csv", sep="")
  raw <- read.csv(filename)[,c("position","count_1","gene")]
  data.frame(position=raw$position, count=raw[,"count_1"],
             gene=raw$gene, pool=pool)
}

reads <- data.frame()
for (pool in pool.names) {
  reads <- rbind(reads, load.pool(pool))
}
READ.CUTOFF <- 15
reads <- reads[reads$count >= READ.CUTOFF,]
positions <- unique(reads$position)
n.positions <- length(positions)
trans <- matrix(F, nrow=n.pools, ncol=n.positions)
row.names(trans) <- pool.names
for (i in 1:n.positions) {
  hits <- reads$pool[reads$position == positions[i]]
  trans[hits,i] <- T
}
row.names(trans) <- NULL

key <- find.closest(pools, trans, status=T)
