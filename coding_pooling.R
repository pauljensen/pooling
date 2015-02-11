G.wagner.24.12.6 <- as.matrix(read.table("wagner_24_14_6.txt"))

expand.G <- function(G) {
  k <- dim(G)[[1]]
  n <- dim(G)[[2]]
  Nk <- 2^k
  u <- matrix(0,nrow=Nk,ncol=k)
  for (i in 0:(Nk-1))
    u[i+1,] <- as.integer(intToBits(i))[1:k]
  (u %*% G) %% 2
}

U <- expand.G(G.wagner.24.12.6)
I <- as.integer(U %*% 2^(0:(dim(U)[2]-1)))
n.ones <- apply(U,1,sum)

count.collisions <- function(ints) {
  collisions <- integer(length(ints))
  for (i in 1:length(ints)) {
    ored <- bitwOr(ints[i], ints)
    collisions[i] <- sum(ored %in% ints)
    print(i)
  }
  collisions
}

collisions <- count.collisions(I)

low <- log(collisions) < 4
U.low <- U[low,]
I.low <- I[low]
n.ones.low <- n.ones[low]
collisions.low <- count.collisions(I.low)


plot(n.ones.low, log(collisions.low))
plot(n.ones, log(collisions))

code.strings <- apply(U.low, 1, function(x) paste(x, sep="", collapse=""))
n.plates <- 33
n.wells <- n.plates * 384
plates <- c(lapply(1:n.plates, function(x) rep(x,384)), recursive=T)
wells <- rep(1:384, n.plates)

codes <- data.frame(plate=plates, well=wells, code=code.strings[1:n.wells])
codes$code <- sample(codes$code, n.wells, replace=F)

# WARNING: overwriting this file without controlling the RNG seed will create
# new pooling codes that differ from the ones used for the experiment.
#write.csv(codes, file="codes.csv",quote=F)


