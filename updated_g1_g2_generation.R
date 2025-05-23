# this function from Koyejo Lab Github https://github.com/koyejo-lab/JointGraphicalLasso
# input is list of matrices same as you would put into `JGL` function, as well as results of `JGL`
# output is AIC for particular estimated graph.
# use this function when tuning lambdas
vBIC <- function(X, est_graph, thr=0.001){
  num <- length(X)
  BIC_acc <- 0.
  for(i in 1:num){
    
    data_num <- dim(X[[i]])[1]
    sample_cov <- cov(X[[i]], X[[i]])
    tr_sum <- sum(diag(sample_cov %*% est_graph$theta[[i]]))
    
    log_det <- determinant(est_graph$theta[[i]], logarithm = TRUE)$modulus[1][1]
    
    E <- sum(sum(abs(est_graph$theta[[i]]) >= thr))
    BIC_acc <- BIC_acc + (tr_sum - log_det) + (log(data_num)*E/data_num)
  }
  return(BIC_acc)
}

# function to make the plots I want below. input data, penalty type and lambdas.
# made function to shorten code later.
lambda_plots <- function(data, penalty = "group", lambda1, lambda2, ...){
  res <- JGL(Y=data,
             penalty=penalty,
             lambda1=lambda1,
             lambda2=lambda2, 
             return.whole.theta=TRUE)
  test1 <- vBIC(data, res, thr=0.0001)
  resgraph <- graph_from_adjacency_matrix(
    -cov2cor(res$theta[[1]]),
    weighted = T,
    mode = "undirected",
    diag = FALSE
  )
  resplot <- plot_jgl(resgraph,
                      multiplier = 3,
                      sub = paste("AIC", round(test1, 2)),
                      ...
  )
}

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

# number of iterations changes how different graphs are?
g2 <- rewire(g1, keeping_degseq(niter = 50))
G_diff_graph <- igraph::union(igraph::difference(g1, g2), 
                              igraph::difference(g2, g1))
gsize(G_diff_graph)
g2_weight <- assignRealDataEdgeWeights(g2, ehist = realDataHist)
g2_weight_boost = boost(g2_weight[[2]])
print(isposdef(as.matrix(g2_weight_boost)))

# get non-zero values and make adjacency matrix
test <- as.matrix(as_adjacency_matrix(g1))

G1 <- test != 0
G1_graph <- graph_from_adjacency_matrix(G1, mode = "undirected")
layout <- layout_nicely(G1_graph)
plot(G1_graph, vertex.frame.color = "white", 
     layout = layout,
     vertex.label = NA, vertex.color = "dodgerblue4", vertex.size = vsize,
     edge.width = 3)



thisSample = sampleNetworkData(N=1000,covMat=solve(g1_weight_boost)) # sample network data function includes standardizing
thisTestSample = sampleNetworkData(N=100,covMat=solve(g1_weight_boost))

sample2 <- sampleNetworkData(N=1000, covMat = solve(g2_weight_boost))

library(pheatmap)
pheatmap(sample2, 
         cluster_cols = T, 
         show_colnames = FALSE, 
         show_rownames = FALSE)

library(JGL)
fgl_results = JGL(Y = list(thisSample, sample2),
                  penalty = "fused",
                  lambda1 = .15,
                  lambda2 = 0.2,
                  return.whole.theta = FALSE)

str(fgl_results) # the theta contains a list of estimated matrices, one for each of the K classes. We will extract the thetas for visualization with igraph. 
print.jgl(fgl_results)

# now tuning lambda

interval_l = 10
lambda.eff <- seq(0.01, 0.4, len = interval_l)
aic_vec1 <- matrix(NA, length(lambda.eff), 1) #for lambda1
aic_vec2 <- matrix(NA, length(lambda.eff), 1) #for lambda2


#search length of lambda1, keeping lambda2 constant and small
start.time = Sys.time()
for(i in 1:length(lambda.eff)){
  fit00 <- JGL(Y=list(thisSample, sample2),penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[1], return.whole.theta=TRUE)
  aic_vec1[i,1] <- vBIC(list(thisSample, sample2), fit00, thr=0.0001)
}

# identify which had min aic
i_idx <- which(aic_vec1 == min(aic_vec1), arr.ind = TRUE)[1,1]
#assign tuned lambda1 based on your search
lam_1 <- lambda.eff[i_idx]

#now search for lambda2, using value of lambda1 learned above
for(j in 1:length(lambda.eff)){
  fit00 <- JGL(Y=list(thisSample, sample2),penalty="group",lambda1=lambda.eff[i_idx],lambda2=lambda.eff[j], return.whole.theta=TRUE)
  aic_vec2[j,1] <- vBIC(list(thisSample, sample2), fit00, thr=0.0001)
}

# identify which lambda2 had min aic
j_idx <- which(aic_vec2 == min(aic_vec2), arr.ind = TRUE)[1,1]
#assign tuned lambda1 based on your search
lam_2 <- lambda.eff[j_idx]

end.time = Sys.time()

end.time-start.time

print(paste("lambda_1", lam_1))
print(paste("lambda_2", lam_2))

fgl_tuned_results <- JGL(Y=list(thisSample, sample2),
                         penalty="fused",
                         lambda1= lam_1,
                         lambda2= lam_2,
                         return.whole.theta=TRUE)

print(fgl_tuned_results)

set.seed(1219)
G <- 4
p <- 51
penalty <- "fused"
lam1 <- 1
lam2 <- 1
v1 <- 1
lam.eff <- lam1 + c(1:10) * 5
v0s <- lam1/lam.eff

start.time = Sys.time()
ssjgl_res <- ssjgl(Y=list(thisSample, sample2),penalty=penalty,lambda0=1, lambda1=lam1,
                   lambda2=lam2, v1 = v1, v0s = v0s, a=1, b=p, doubly=TRUE, normalize=FALSE)
end.time = Sys.time()
end.time-start.time

# extract all estimated covariance matrices from result
ssjgl_covar_matrices <- ssjgl_res$thetalist[[length(v0s)]]
ss_res1 <- (abs(ssjgl_covar_matrices[[1]])>0.001)*1
diag(ss_res1) = 0
ss_perf1 <- evaluatePerformance(G1, ss_res1)

ss_res_df <- data.frame(
  tpr = ss_perf1[[1]]/(ss_perf1[[1]] + ss_perf1[[4]]),
  fpr = ss_perf1[[3]]/(ss_perf1[[3]] + ss_perf1[[2]]),
  auc = auc(roc(G1~ss_res1))
  
)

# get 5% percentile of absolute value of empirical dist
thresh <- quantile(abs(realDataHist$mids), .001)

curr_res1 <-  (abs(fgl_tuned_results$theta[[1]]) > thresh)*1

diag(curr_res1) = 0


curr_perf1 <- evaluatePerformance(G1, curr_res1)

curr_res2 <-  (abs(fgl_tuned_results$theta[[2]]) > thresh)*1

diag(curr_res2) = 0

library(pROC)
roc_res <- roc(G1~curr_res1)
auc(roc_res)


curr_perf2 <- evaluatePerformance(G2, curr_res2)


res_df <- data.frame(
  tpr = curr_perf1[[1]]/(curr_perf1[[1]] + curr_perf1[[4]]),
  fpr = curr_perf1[[3]]/(curr_perf1[[3]] + curr_perf1[[2]]),
  auc = auc(roc_res),
  
)
