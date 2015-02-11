
for (thresh in 0:4) {
  pos <- positions[key$distance <= thresh]
  genes <- unique(reads$gene[reads$position %in% pos])
  print(c(length(pos), length(genes)))
}


# export mapping
thresh <- 3
loc <- key$pool.index[key$distance <= thresh]
pos <- positions[key$distance <= thresh]
lookup <- cbind(codes[loc,], key[key$distance <= thresh,],
                data.frame(position=pos,gene=""))
lookup$gene <- as.character(lookup$gene)
for (i in 1:nrow(lookup)) {
  genes <- reads$gene[reads$position == lookup$position[i]]
  lookup$gene[i] <- as.character(genes)[1]
}

write.csv(lookup, "lookup_thresh3.csv")
