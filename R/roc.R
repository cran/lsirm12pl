
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

  auc = round(roc_score$auc,4)
  legend_text <- paste0("AUC = ", round(auc, 2))
  rocp = ggroc(roc_score, colour = 'steelblue', size = 2)+
    theme(axis.text.x = element_text(face="bold",size=15),
          axis.text.y = element_text( face="bold",size=15),
          axis.title=element_text(size=15, face='bold'),
          plot.margin = margin(1,1,1.5,1.2,"cm"))+
    annotate("text", x = 0.3, y = 0.05, label = legend_text,
             fontface=2, size=4) +
    xlab("Specificity")+
    ylab("Sensitivity")
  return(rocp)
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

  auc = round(roc_score$auc,4)
  legend_text <- paste0("AUC = ", round(auc, 2))
  rocp = ggroc(roc_score, colour = 'steelblue', size = 2)+
    theme(axis.text.x = element_text(face="bold",size=15),
          axis.text.y = element_text( face="bold",size=15),
          axis.title=element_text(size=15, face='bold'),
          plot.margin = margin(1,1,1.5,1.2,"cm"))+
    annotate("text", x = 0.3, y = 0.05, label = legend_text,
             fontface=2, size=4)+
    xlab("Specificity")+
    ylab("Sensitivity")
  return(rocp)
}


