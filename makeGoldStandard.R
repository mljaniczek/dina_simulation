library(ggplot2)
library(huge)
library(igraph)
library(MASS)
library(moments)
library(pracma)

# read in reference weights
realDataHist = read.csv("mxDist.csv")

# sample network data
sampleNetworkData = function(N, covMat)
{
  sample = mvrnorm(N, mu = rep(0,nrow(covMat)), covMat)
  return(sample)
}


#### Assign Edge Weights ####

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

#### Make Gold Standard Network Structures #### 

makeGoldStandardNets = function(P, reference_file)
{
  PminusOne = P-1
  
  # make random graph
  erHigh = sample_gnp(P,0.20)
  erHighDensity = length(E(erHigh))/(P*PminusOne/2)
  erLow = sample_gnp(P,0.06)
  erLowDensity = length(E(erLow))/(P*PminusOne/2)
  
  print(erHighDensity)
  print(erLowDensity)
  
  # make scale-free graph
  sfHigh = sample_pa(n=P, power=1, m=5, directed=F)
  sfHighDensity = length(E(sfHigh))/(P*PminusOne/2)
  sfLow = sample_pa(n=P, power=1, m=2, directed=F)
  sfLowDensity = length(E(sfLow))/(P*PminusOne/2)
  
  print(sfHighDensity)
  print(sfLowDensity)
  
  erHighRealData = assignRealDataEdgeWeights(erHigh,ehist=reference_file)
  erLowRealData = assignRealDataEdgeWeights(erLow,ehist=reference_file)
  sfHighRealData = assignRealDataEdgeWeights(sfHigh,ehist=reference_file)
  sfLowRealData = assignRealDataEdgeWeights(sfLow,ehist=reference_file)
  
  erLowPrecRealDataBoost = boost(erLowRealData[[2]])
  erHighPrecRealDataBoost = boost(erHighRealData[[2]])
  sfLowPrecRealDataBoost = boost(sfLowRealData[[2]])
  sfHighPrecRealDataBoost = boost(sfHighRealData[[2]])
  
  adjMatListRealData = c(erLowPrecRealDataBoost,
                         erHighPrecRealDataBoost,
                         sfLowPrecRealDataBoost,
                         sfHighPrecRealDataBoost)
  
  
  for(mat in adjMatListRealData)
  {
    print(isposdef(as.matrix(mat)))
  }
  
  return(adjMatListRealData)
}

test <- makeGoldStandardNets(50, realDataHist)

test2 <- as.matrix(test[[1]])
t1 <- test2 != 0

test3 <- rewire(graph_from_adjacency_matrix(t1, mode = "undirected"), keeping_degseq(niter = 30))


graphLabels = c("Random - Low Density",
                "Random - High Density",
                "Scale-Free - Low Density",
                "Scale-Free - High Density")

pdf("graphLayouts.pdf",width=8,height=16)
par(mfrow=c(4,2),mar=c(1,1,1,1))
for(m in 1:4)
{
  thisGraph = graph_from_adjacency_matrix(-cov2cor(as.matrix(test[[m]])),
                                          weighted=T,
                                          mode="undirected",
                                          diag=F)
  dummyGraph = thisGraph
  E(dummyGraph)$weight = NA #ifelse(E(thisGraph)$weight > 0, 1, 0)
  myLayout = layout_with_fr(dummyGraph)
  plot(thisGraph, vertex.label=NA,
       edge.width = ifelse(abs(E(thisGraph)$weight)>0,1.5,0),
       vertex.size = (degree(thisGraph)+2),
       vertex.color = rainbow(50)[degree(thisGraph)],
       layout=myLayout,
       main = graphLabels[m])
}
 dev.off()