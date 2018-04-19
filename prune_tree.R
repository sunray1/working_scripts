intree <- read.nexus("test_tree.tre")
iterations <- 1


for (i in 1:iterations){
  dists <- distRoot(intree)
  cutoff <- mean(dists) + sd(dists)*2.5
  long_dists <- dists[dists > cutoff]
  intree <- drop.tip(intree, names(long_dists))
  write.nexus(intree, file = paste("test_tree_",i,".tre", sep=""))
}
