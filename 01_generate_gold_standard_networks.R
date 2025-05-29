library(pracma) # for isposdef


realDataHist = read.csv("mxDist.csv")
source("00_helper_functions.R")

library(MatrixCorrelation)

relativeFrobNormAfter = function(estimatedMatrix, trueMatrix)
{
  Delta = as.matrix(estimatedMatrix - trueMatrix)
  return(norm(Delta,"F")/norm(as.matrix(trueMatrix),"F"))
}

# input g1 adjacency matrix generated from igraph. can be random, scale-free etc.
# perc is percent to change in rewire

gen_gold <- function(g1, perc, empirical_density){

  
  g2 = rewire(g1, keeping_degseq(niter = gsize(g1)*perc))
  
  # before weighting 
  g1_unif <- g1_emp <- g1
  g2_unif <- g2_emp <- g2
  
  # assign edge weights using emprical distribution
  weights = sample(empirical_density$mids, replace=T,size=length(E(g1)), prob=empirical_density$density)  
  weights1 = weights
  weights2 = weights
  
  # use this if you want to randomly make some of the weights 0 in each
  #new_weight_index1 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.05, .95)) 
  #new_weight_index2 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.05, .95)) 
  #weights1[new_weight_index1] <- 0
  #weights2[new_weight_index2] <- 0
  
  E(g1_emp)$weights <- weights1
  E(g2_emp)$weights <- weights2
  
  
  # assign edge weights sampling from uniform distribution
  size <- length(E(g1))
  a <- 0.01
  b <-  0.75
  samp_right <- runif(ceiling(size / 2), min = a, max = b)
  samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
  weights_unif <- sample(c(samp_left, samp_right), size)
  weights_unif1 <- weights_unif2 <- weights_unif
  #weights_unif1[new_weight_index1] <- 0
  #weights_unif2[new_weight_index2] <- 0

  E(g1_unif)$weights <- weights_unif1
  E(g2_unif)$weights <- weights_unif2
  
  precMat1_emp = boost(-1*as_adjacency_matrix(g1_emp, attr = c("weights"), type="both") + diag(length(g1_emp)))
  precMat2_emp = boost(-1*as_adjacency_matrix(g2_emp, attr = c("weights"), type="both") + diag(length(g2_emp)))
  precMat1_unif = boost(-1*as_adjacency_matrix(g1_unif, attr = c("weights"), type="both") + diag(length(g1_unif)))
  precMat2_unif = boost(-1*as_adjacency_matrix(g2_unif, attr = c("weights"), type="both") + diag(length(g2_unif)))

  dat = list(
    g = list("g1" = g1, "g2" = g2),
    graphs = list("g1_emp" = g1_emp, 
                  "g2_emp" = g2_emp, 
                  "g1_unif" = g1_unif, 
                  "g2_unif" = g2_unif),
    omega = list("o1_emp" = precMat1_emp, 
                 "o2_emp" = precMat2_emp,
                 "o1_unif" = precMat1_unif,
                 "o2_unif" = precMat2_unif),
    sigma = list("s1_emp" = solve(precMat1_emp),
                 "s2_emp" = solve(precMat2_emp),
                 "s1_unif" = solve(precMat1_unif),
                 "s2_unif" = solve(precMat2_unif))
    )
  
  return(dat)
}

#test <- gen_gold(g1, .1, realDataHist)


# function to generate data from gold standard graphs and return list of either empirical or uniform weighted 
gen_dat <- function(gold, n){
  dat1_emp = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o1_emp)), Sigma = gold$sigma$s1_emp))
  dat2_emp = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o2_emp)), Sigma = gold$sigma$s2_emp))
  dat1_unif = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o1_unif)), Sigma = gold$sigma$s1_unif))
  dat2_unif = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o2_unif)), Sigma = gold$sigma$s2_unif))
  
  dat = list(
    "emp_dat" = list("dat1" = dat1_emp,
                     "dat2" = dat2_emp),
    "unif_dat" = list("dat1" = dat1_unif,
                      "dat2" = dat2_unif)
  )
  
  return(dat)
}
# 
# test_dat <- gen_dat(test, 1000)
# 
# # function to tune lambdas in JGL
# # now tuning lambda
# 
# 
# # get 5% percentile of absolute value of empirical dist
# thresh <- 0.005 #quantile(abs(weights), .025)

tune_jgl <- function(dat_list, threshold, truth, true_omega1, true_omega2, setting = "setting"){
  g1 = truth$g$g1
  g2 = truth$g$g2
  X = dat_list
  interval_l = 10
  lambda.eff <- seq(0.01, 0.4, len = interval_l)
  aic_vec1 <- matrix(NA, length(lambda.eff), 1) #for lambda1
  aic_vec2 <- matrix(NA, length(lambda.eff), 1) #for lambda2

    #search length of lambda1, keeping lambda2 constant and small
    start.time = Sys.time()
    for(i in 1:length(lambda.eff)){
      fit00 <- JGL(Y=X,penalty="fused",lambda1=lambda.eff[i],lambda2=lambda.eff[1], return.whole.theta=TRUE)
      aic_vec1[i,1] <- vBIC(X, fit00, thr=0.0001)
    }
    
    # identify which had min aic
    i_idx <- which(aic_vec1 == min(aic_vec1), arr.ind = TRUE)[1,1]
    #assign tuned lambda1 based on your search
    lam_1 <- lambda.eff[i_idx]
    
    #now search for lambda2, using value of lambda1 learned above
    for(j in 1:length(lambda.eff)){
      fit00 <- JGL(Y=X,penalty="group",lambda1=lambda.eff[i_idx],lambda2=lambda.eff[j], return.whole.theta=TRUE)
      aic_vec2[j,1] <- vBIC(X, fit00, thr=0.0001)
    }
    
    # identify which lambda2 had min aic
    j_idx <- which(aic_vec2 == min(aic_vec2), arr.ind = TRUE)[1,1]
    #assign tuned lambda1 based on your search
    lam_2 <- lambda.eff[j_idx]
    
    end.time = Sys.time()
    
    print(end.time-start.time)

    # 6.33 minutes for 500 features, 1000 per arm
    # 48 seconds for 200 features, 1000 per arm
    
    res_fgl <- JGL(Y=X,
                   penalty="fused",
                   lambda1= lam_1,
                   lambda2= lam_2,
                   return.whole.theta=TRUE)

#save(res_fgl, file = "fgl_res_n1000p200.rdata")

      
      
      curr_res1 <-  (abs(res_fgl$theta[[1]]) > thresh)*1
      
      diag(curr_res1) = 0
      
      
      curr_perf1 <- evaluatePerformance(g1, curr_res1)
      
      curr_res2 <-  (abs(res_fgl$theta[[2]]) > thresh)*1
      
      diag(curr_res2) = 0
      
      
      curr_perf2 <- evaluatePerformance(g2, curr_res2)
      
      
      res_df <- data.frame(
        setting = setting,
        tpr1 = curr_perf1[[1]]/(curr_perf1[[1]] + curr_perf1[[4]]),
        fpr1 = curr_perf1[[3]]/(curr_perf1[[3]] + curr_perf1[[2]]),
        auc1 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(g1)) ~ curr_res1))),
        frob1 = relativeFrobNormAfter(res_fgl$theta[[1]], true_omega1),
        tpr2 = curr_perf2[[1]]/(curr_perf2[[1]] + curr_perf2[[4]]),
        fpr2 = curr_perf2[[3]]/(curr_perf2[[3]] + curr_perf2[[2]]),
        auc2 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(g2)) ~ curr_res2))),
        frob2 = relativeFrobNormAfter(res_fgl$theta[[2]], true_omega2)
      )
      
      return(res_df)
}






# delta <- as.matrix(test$omega$o1_emp)
# relativeFrobNormAfter(res_fgl$theta[[1]], test$omega$o1_emp)
# 
# test_res_emp <- tune_jgl(test_dat$emp_dat, thresh, test, test$omega$o1_emp, test$omega$o2_emp)
# test_res_unif <- tune_jgl(test_dat$unif_dat, thresh, test, test$omega$o1_unif, test$omega$o2_unif)
# #test_dat <- gen_dat(test, 1000)
# 

# now function to run ssjgl and get results

run_ssjgl <- function(dat_list, threshold, truth, true_omega1, true_omega2, setting = "setting"){
  X = dat_list
  g1 = truth$g$g1
  g2 = truth$g$g2
  penalty <- "fused"
  lam1 <- 1
  lam2 <- 1
  v1 <- 1
  lam.eff <- lam1 + c(1:4) * 5
  v0s <- lam1/lam.eff
  
  start.time = Sys.time()
  ssjgl_res <- ssjgl(Y=X,penalty=penalty,lambda0=1, lambda1=lam1,
                     lambda2=lam2, v1 = v1, v0s = v0s, a=1, b=p, doubly=TRUE, normalize=FALSE)
  end.time = Sys.time()
  print(end.time-start.time)
  
  # extract all estimated inverse covariance matrices from result
  ssjgl_covar_matrices <- ssjgl_res$thetalist[[length(v0s)]]
  ss_res1 <- (abs(ssjgl_covar_matrices[[1]])>thresh)*1
  diag(ss_res1) = 0
  ss_perf1 <- evaluatePerformance(g1, ss_res1)
  ss_res2 <- (abs(ssjgl_covar_matrices[[2]])>thresh)*1
  diag(ss_res2) = 0
  ss_perf2 <- evaluatePerformance(g2, ss_res2)
  
  ss_res_df <- data.frame(
    setting = setting,
    tpr1 = ss_perf1[[1]]/(ss_perf1[[1]] + ss_perf1[[4]]),
    fpr1 = ss_perf1[[3]]/(ss_perf1[[3]] + ss_perf1[[2]]),
    auc1 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(g1))~ss_res1))),
    frob1 = relativeFrobNormAfter(ssjgl_covar_matrices[[1]], true_omega1),
              tpr2 = ss_perf2[[1]]/(ss_perf2[[1]] + ss_perf2[[4]]),
              fpr2 = ss_perf2[[3]]/(ss_perf1[[3]] + ss_perf2[[2]]),
              auc2 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(g1))~ss_res2))),
    frob2 = relativeFrobNormAfter(ssjgl_covar_matrices[[2]], true_omega2)
  )
  
  return(ss_res_df)  
  
}

# relativeFrobNormAfter(ssjgl_covar_matrices[[2]], test$omega$o2_emp)
# 
# test_ss_res_emp <- run_ssjgl(test_dat$emp_dat, thresh, test, test$omega$o1_emp, test$omega$o2_emp)
# test_ss_res_unif <- run_ssjgl(test_dat$unif_dat, thresh, test,test$omega$o1_emp, test$omega$o2_emp)
