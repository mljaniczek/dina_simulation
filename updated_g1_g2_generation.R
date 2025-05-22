# assign real edge weights using reference data

assignRealDataEdgeWeights = function(graph,ehist)
{
  P = length(V(graph))
  weights = sample(ehist$mids, replace=T,size=length(E(graph)), prob=ehist$density)  
  E(graph)$weights = weights
  precMat = -1*as_adjacency_matrix(graph, attr = c("weights"), type="both") + diag(P)
  return(list(graph, precMat))
}


#### Adjust Covariance Matrices to be Positive Definite ####

boost = function(myMatrix)
{
  minEig = min(eigen(myMatrix)$values)
  if(minEig < 0)
  {
    print("Boosting!")
    return(myMatrix - minEig*1.01*diag(ncol(myMatrix)))
  }
  else
    return(myMatrix)
}

realDataHist = read.csv("mxDist.csv")

# make random graph
P = 100
dens = 0.20

g1 = sample_gnp(P,dens)
g1Density = length(E(g1))/(P*(P-1)/2)

g1_weight = assignRealDataEdgeWeights(g1,ehist=realDataHist)

g1_weight_boost = boost(g1_weight[[2]])

print(isposdef(as.matrix(g1_weight_boost)))


g2 <- rewire(g1, keeping_degseq(niter = 3))
g2_weight <- assignRealDataEdgeWeights(g2, ehist = realDataHist)
g2_weight_boost = boost(g2_weight[[2]])
print(isposdef(as.matrix(g2_weight_boost)))


ggcorrplot(as.matrix(g1_weight_boost))
ggcorrplot(as.matrix(g2_weight_boost))



