
roc_1pl = function(lsirm_result, ...){
  data = lsirm_result$data
  beta = lsirm_result$beta_estimate
  theta = lsirm_result$theta_estimate
  gamma = lsirm_result$gamma_estimate
  z = lsirm_result$z_estimate
  w = lsirm_result$w_estimate

  p = ncol(data)
  n = nrow(data)
  M1 = rep(1, n) %*% t(beta)
  M2 = theta %*% t(rep(1, p))
  M3 = diag(z %*% t(z)) %*% t(rep(1, p))
  M4 = rep(1, n) %*% t(diag(w %*% t(w)))
  prob = (M1 + M2) - gamma*sqrt(M3 + M4 -2*(z %*% t(w)))
  roc_score = pROC::roc(c(data.matrix(data)),c(prob),
                        levels = c(0,1), direction = "<")
  # roc_score = pROC::roc(c(data.matrix(data)),c(prob), ...)
  plot(roc_score, main ="ROC Curve", print.auc = T)
  return(roc_score$auc)
}

roc_2pl = function(lsirm_result){
  data = lsirm_result$data
  beta = lsirm_result$beta_estimate
  theta = lsirm_result$theta_estimate
  alpha = lsirm_result$alpha_estimate
  gamma = lsirm_result$gamma_estimate
  z = lsirm_result$z_estimate
  w = lsirm_result$w_estimate

  p = ncol(data)
  n = nrow(data)
  M1 = rep(1, n) %*% t(beta)
  M2 = theta %*% t(alpha)
  M3 = diag(z %*% t(z)) %*% t(rep(1, p))
  M4 = rep(1, n) %*% t(diag(w %*% t(w)))
  prob = (M1 + M2) - gamma*sqrt(M3 + M4 -2*(z %*% t(w)))

  roc_score = pROC::roc(c(data.matrix(data)),c(prob),
                        levels = c(0,1), direction = "<")
  plot(roc_score, main ="ROC curve", print.auc = T)
  return(roc_score$auc)
}


