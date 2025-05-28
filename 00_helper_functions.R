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
  minEig = min(as.numeric(eigen(myMatrix)$values))
  if(minEig < 0)
  {
    print("Boosting!")
    return(myMatrix - minEig*1.01*diag(ncol(myMatrix)))
  }
  else
    return(myMatrix)
}
