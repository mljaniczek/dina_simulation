library(pracma) # for isposdef


realDataHist = read.csv("mxDist.csv")
source("00_helper_functions.R")

# make random graph
p = 100
dens2 = 0.06
perc = .10

gen_gold_nets <- function(p, dens1, dens2, weights, perc){
  
  # G1
  # random graph
  erhigh1 = sample_gnp(p, dens1)

  erlow1 = sample_gnp(p, dens2)
  erhigh1_dens = length(E(erhigh1))/p*(p-1)/2
  erhigh1_dens = length(E(erhigh2))/p*(p-1)/2
  
  # make scale-free graph
  sfHigh1 = sample_pa(n=p, power=1, m=5, directed=F)
  sfHighDensity = length(E(sfHigh1))/(p*(p-1)/2)
  sfLow1 = sample_pa(n=p, power=1, m=2, directed=F)
  sfLowDensity = length(E(sfLow1))/(p*(p-1)/2)

  # G2
  erhigh2 = rewire(erhigh1, keeping_degseq(niter = gsize(erhigh1)*perc))
  erlow2 = rewire(erlow1, keeping_degseq(niter = gsize(erlow1)*perc))
  
  G_diff_graph1 <- igraph::union(igraph::difference(erhigh1, erhigh2), 
                                igraph::difference(erhigh2, erhigh1))
  G_diff_graph2 <- igraph::union(igraph::difference(erlow1, erlow2), 
                                 igraph::difference(erlow2, erlow1))
  gsize(G_diff_graph1)
  gsize(G_diff_graph2)
  
  
  # assign edge weights
  weights = sample(realDataHist$mids, replace=T,size=length(E(erhigh1)), prob=realDataHist$density)  
  weights1 = weights
  weights2 = weights
  new_weight_index1 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.1, .9)) 
  new_weight_index2 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.1, .9)) 
  weights1[new_weight_index1] <- 0
  weights2[new_weight_index2] <- 0
  
  E(erhigh1)$weights <- weights1
  E(erhigh2)$weights <- weights2
  
  size <- length(E())
  a <- 0.2
  b <-  0.9
  samp_right <- runif(ceiling(size / 2), min = a, max = b)
  samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
  weights_unif <- sample(c(samp_left, samp_right), size)
  
  
  
  return(dat)
}


g1 = sample_gnp(p,dens)

g1_sf = sample_pa(n=p, power=1, m=5, directed=F)

g1Density = length(E(g1))/(p*(p-1)/2)
g1_sf_dens = length(E(g1_sf))/(p*(p-1)/2)

g1_weight = assignRealDataEdgeWeights(g1,ehist=realDataHist)

g1_weight_boost = boost(g1_weight[[2]])

print(isposdef(as.matrix(g1_weight_boost)))

# number of iterations changes how different graphs are?
g2 <- rewire(g1, keeping_degseq(niter = 50))
G_diff_graph <- igraph::union(igraph::difference(g1, g2), 
                              igraph::difference(g2, g1))
gsize(G_diff_graph)
g2_weight <- assignRealDataEdgeWeights(g2, ehist = realDataHist)
g2_weight_boost = boost(g2_weight[[2]])
print(isposdef(as.matrix(g2_weight_boost)))