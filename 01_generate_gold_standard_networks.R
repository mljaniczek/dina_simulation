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

gen_gold <- function(g1, perc, empirical_density, wt_method){

  
  g2 = rewire(g1, keeping_degseq(niter = gsize(g1)*perc))
  g_diff_graph <- igraph::union(igraph::difference(g1, g2), 
                                igraph::difference(g2, g1))
  
  g1_graph <- g1
  g2_graph <- g2
  
  if(wt_method == "empirical"){
    # assign edge weights using emprical distribution
    weights = sample(empirical_density$mids, replace=T,size=length(E(g1)), prob=empirical_density$density)  
    weights1 = weights
    weights2 = weights
    # use this if you want to randomly make some of the weights 0 in each
    #new_weight_index1 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.05, .95)) 
    #new_weight_index2 <- sample(c(TRUE, FALSE), replace = T, size = length(weights1), prob = c(.05, .95)) 
    #weights1[new_weight_index1] <- 0
    #weights2[new_weight_index2] <- 0
    E(g1_graph)$weights <- weights1
    E(g2_graph)$weights <- weights2
  } else{
    # assign edge weights sampling from uniform distribution
    size <- length(E(g1))
    a <- 0.035
    b <-  0.75
    samp_right <- runif(ceiling(size / 2), min = a, max = b)
    samp_left <- runif(ceiling(size / 2), min = -b, max = -a)
    weights_unif <- sample(c(samp_left, samp_right), size)
    weights_unif1 <- weights_unif2 <- weights_unif
    #weights_unif1[new_weight_index1] <- 0
    #weights_unif2[new_weight_index2] <- 0
    
    E(g1_graph)$weights <- weights_unif1
    E(g2_graph)$weights <- weights_unif2
  }
  precMat1 = boost(-1*as_adjacency_matrix(g1_graph, attr = c("weights"), type="both") + diag(length(g1_graph)))
  precMat2 = boost(-1*as_adjacency_matrix(g2_graph, attr = c("weights"), type="both") + diag(length(g2_graph)))
  
  diff_true = (abs(precMat1-precMat2) > 0.1)*1

  dat = list(
    g = list("g1" = g1, "g2" = g2),
    graphs = list("g1_graph" = g1_graph, 
                  "g2_graph" = g2_graph),
    omega = list("o1" = precMat1, 
                 "o2" = precMat2),
    sigma = list("s1" = solve(precMat1),
                 "s2" = solve(precMat2)),
    g_diff_graph = diff_true,
    g_diff_graph2 = g_diff_graph
    )
  
  return(dat)
}




# function to generate data from gold standard graphs and return list of either empirical or uniform weighted 
gen_dat <- function(gold, n){
  dat1 = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o1)), Sigma = gold$sigma$s1))
  dat2 = scale(mvrnorm(n = n, mu = rep(0,nrow(gold$omega$o2)), Sigma = gold$sigma$s2))
  dat = list("dat1" = dat1,
              "dat2" = dat2)
  
  return(dat)
}
# 

# 
# # function to tune lambdas in JGL
# # now tuning lambda
# 
# 
# # get 5% percentile of absolute value of empirical dist
# thresh <- 0.005 #quantile(abs(weights), .025)

tune_jgl <- function(dat_list){
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
    


      return(list("lam_1" = lam_1, "lam_2" = lam_2))
}


# calculate error metrics 
calc_metrics <- function(thetas, truth, threshold, setting, diff_thresh){
  
  curr_res1 <-  (abs(thetas[[1]]) > threshold)*1
  
  diag(curr_res1) = 0
  
  curr_perf1 <- suppressMessages(evaluatePerformance(truth$g$g1, curr_res1))
  
  curr_res2 <-  (abs(thetas[[2]]) > threshold)*1
  
  diag(curr_res2) = 0
  
  curr_perf2 <- suppressMessages(evaluatePerformance(truth$g$g2, curr_res2))
  
 # diff_res <- abs(curr_res1-curr_res2)
 # diag(diff_res) <- 0
  
  #diff_perf <- suppressMessages(evaluatePerformance(as_adjacency_matrix(truth$g_diff_graph2), diff_res))
  
  diff_res2 <- (abs(thetas[[1]] - thetas[[2]]) > diff_thresh)*1
  diag(diff_res2) = 0
  diff_perf2 <- suppressMessages(evaluatePerformance(truth$g_diff_graph, diff_res2))
  
  
  res_df <- data.frame(
    setting = setting,
    tpr1 = curr_perf1[[1]]/(curr_perf1[[1]] + curr_perf1[[4]]),
    fpr1 = curr_perf1[[3]]/(curr_perf1[[3]] + curr_perf1[[2]]),
    auc1 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(truth$g$g1)) ~ curr_res1))),
    #frob1 = relativeFrobNormAfter(thetas[[1]], truth$omega$o1),
    tpr2 = curr_perf2[[1]]/(curr_perf2[[1]] + curr_perf2[[4]]),
    fpr2 = curr_perf2[[3]]/(curr_perf2[[3]] + curr_perf2[[2]]),
    auc2 = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(truth$g$g2)) ~ curr_res2))),
    #frob2 = relativeFrobNormAfter(thetas[[2]], truth$omega$o2)
    # tpr_diff = diff_perf[[1]]/(diff_perf[[1]] + diff_perf[[4]]),
    # fpr_diff = diff_perf[[3]]/(diff_perf[[3]] + diff_perf[[2]]),
   #  auc_diff = as.numeric(auc(roc(as.matrix(as_adjacency_matrix(truth$g_diff_graph2)) ~ diff_res))),
    tpr_diff2 = diff_perf2[[1]]/(diff_perf2[[1]] + diff_perf2[[4]]),
    fpr_diff2 = diff_perf2[[3]]/(diff_perf2[[3]] + diff_perf2[[2]]),
    auc_diff2 = as.numeric(auc(roc(as.matrix(truth$g_diff_graph) ~ diff_res2)))
  )  %>%
    mutate_if(is.numeric, round, 3)
  return(res_df)
}

get_results <- function(n, p, gold_list_uniform, gold_list_empirical, v0s, t_num){
  
  #generate data
  
  dat_list_uniform <- map(gold_list_uniform, ~gen_dat(.x, n))
  dat_list_empirical <- map(gold_list_empirical, ~gen_dat(.x, n))
  
  # just run this tuning on one dataset
 # lam_list_emp <- tune_jgl(dat_list_empirical$random_low)
  lam_list_emp <- map(dat_list_empirical, ~tune_jgl(.x))
  lam_list_unif <- map(dat_list_uniform, ~tune_jgl(.x))
  
  jgl_res_unif <- map2(dat_list_uniform, lam_list_unif, ~JGL(Y=.x,
                                                             penalty="fused",
                                                             lambda1= .y$lam_1,
                                                             lambda2= .y$lam_2,
                                                             return.whole.theta=TRUE))
  
  jgl_res_emp <- map2(dat_list_empirical, lam_list_emp, ~JGL(Y=.x,
                                                             penalty="fused",
                                                             lambda1= .y$lam_1,
                                                             lambda2= .y$lam_2,
                                                             return.whole.theta=TRUE))
  
  tic()
  ssjgl_res_emp <- map(dat_list_empirical, ~ ssjgl(Y=.x,penalty="fused",lambda0=1, lambda1=1,
                                                   lambda2=1, v1 = 1, v0s = v0s, a=1, b=1, doubly=TRUE, normalize=FALSE))
  toc()
  
  tic()
  ssjgl_res_unif <- map(dat_list_uniform, ~ ssjgl(Y=.x,penalty="fused",lambda0=1, lambda1=1,
                                                  lambda2=1, v1 = 1, v0s = v0s, a=1, b=1, doubly=TRUE, normalize=FALSE))
  toc()
  
  
  
  jgl_res_err_unif <- map2(jgl_res_unif,gold_list_uniform, ~calc_metrics(.x$theta, .y, 0.01, "uniform", 0.03)) %>%
    bind_rows(.id = "data") %>%
    mutate(method = "jgl")
  
  ssjgl_res_err_unif <- map2(ssjgl_res_unif, gold_list_uniform, ~calc_metrics(.x$thetalist[[t_num]], .y, .01, "uniform", 0.03))%>%
    bind_rows(.id = "data") %>%
    mutate(method = "sjgl")
  
  jgl_res_err_emp <- map2(jgl_res_emp,gold_list_empirical, ~calc_metrics(.x$theta, .y, 0.01, "emp", 0.03))%>%
    bind_rows(.id = "data") %>%
    mutate(method = "jgl")
  
  ssjgl_res_err_emp <- map2(ssjgl_res_emp, gold_list_empirical, ~calc_metrics(.x$thetalist[[t_num]], .y, .01, "emp", 0.03))%>%
    bind_rows(.id = "data") %>%
    mutate(method = "sjgl")
  
  res_master <- bind_rows(jgl_res_err_unif, ssjgl_res_err_unif, jgl_res_err_emp, ssjgl_res_err_emp) %>%
    mutate(p = p,
           n = n)
  
  return(res_master)

}

# 
# nab_error_mets <- function(res, t_num, gold_list_uniform, gold_list_empirical){
#   
# }
 # delta <- as.matrix(test$omega$o1_emp)
# relativeFrobNormAfter(res_fgl$theta[[1]], test$omega$o1_emp)
# 
# test_res_emp <- tune_jgl(test_dat$emp_dat, thresh, test, test$omega$o1_emp, test$omega$o2_emp)
# test_res_unif <- tune_jgl(test_dat$unif_dat, thresh, test, test$omega$o1_unif, test$omega$o2_unif)
# #test_dat <- gen_dat(test, 1000)
# 

# now function to run ssjgl and get results

tune_ssjgl <- function(dat_list){
  X = dat_list
  penalty <- "fused"
  lam1 <- 1
  lam2 <- 1
  v1 <- 1
  lam.eff <- lam1 + c(1:10) * 5
  v0s <- lam1/lam.eff
  
  start.time = Sys.time()
  ssjgl_res <- ssjgl(Y=X,penalty=penalty,lambda0=1, lambda1=lam1,
                     lambda2=lam2, v1 = v1, v0s = v0s, a=1, b=1, doubly=TRUE, normalize=FALSE)
  end.time = Sys.time()
  print(end.time-start.time)
  
  # extract all estimated inverse covariance matrices from result
  ssjgl_covar_matrices <- ssjgl_res$thetalist[[length(v0s)]]
  
  return(ssjgl_res)  
  
}


# K <- 3
# p <- 50
# n <- 20
# data <- generateData_rewire(K = K, p = p, n = n, ncores = 1, verbose = FALSE)
# G_common_true <- data$CommonGraph
# X <- data$Data
# res <- jewel(X, lambda1 = 0.25)
# G_common_est <- res$CommonG
# evaluatePerformance(G = G_common_true, G_hat = G_common_est)

# test 
# 
# test_gold <- gen_gold(g1_random_high, .1, new_realdat, wt_method = "uniform")
# test_dat <- gen_dat(test_gold, 1000)
# test <- tune_jgl(test_dat)
# test_ss <- tune_ssjgl(test_dat)
# 
# res_fgl <- JGL(Y=test_dat,
#                penalty="fused",
#                lambda1= test$lam_1,
#                lambda2= test$lam_2,
#                return.whole.theta=TRUE)
# 
# test_res_err <- calc_metrics(res_fgl$theta, test_gold, 0.01, "random_high_uniform", .1)
# 
# 
# test_res_err_ss3 <- calc_metrics(test_ss$thetalist[[4]], test_gold, 0.01, "random_high_uniform", .06)


# relativeFrobNormAfter(ssjgl_covar_matrices[[2]], test$omega$o2_emp)
# 
# test_ss_res_emp <- run_ssjgl(test_dat$emp_dat, thresh, test, test$omega$o1_emp, test$omega$o2_emp)
# test_ss_res_unif <- run_ssjgl(test_dat$unif_dat, thresh, test,test$omega$o1_emp, test$omega$o2_emp)
