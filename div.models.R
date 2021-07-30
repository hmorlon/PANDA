# This function mimics div.models from functions.for.shift.estimates.R but with RPANDA functions
# backbone options are available.

div.models <- function(phylo, tot_time, f,
                       backbone = F, spec_times = NULL, branch_times = NULL,
                       models = c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR"),
                       cond, verbose = T, n.max = NULL, rate.max = NULL){
  
  results <- as.data.frame(matrix(NA,length(models),8))
  colnames(results)<-c("Models","Parameters","logL","AICc","Lambda","Alpha","Mu","Beta")
  res_l <- 1
  
  list_param<- list(c(0.1, 0.01, 0.02, 0.01))
  list_param[2] <- list(c(0.1, 0.01, 0.02, 0.01))
  list_param[3] <- list(c(0.05, 0.01, 0.01, 0.01))
  list_param[4] <- list(c(0.1, 0.01, 0.005, 0.01))
  list_param[5] <- list(c(0.01, 0.001, 0.001, 0.001))
  list_param[6] <- list(c(0.3, 0.01, 0.05, 0.01))
  
  names(list_param) <- c("improved", "base1", "base2", "base3", "base4", "base5")
  
  if("BCST" %in% models){
    
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        suppressWarnings(test <- try(
          treei_BCST <- fit_bd_backbone(phylo, tot_time, f=f,
                                        backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                        f.lamb = function(x,y){y}, f.mu = function(x,y){0},
                                        lamb_par = list_param[[p]][1], mu_par = NULL,
                                        cst.lamb=T, cst.mu=T, expo.lamb=F, expo.mu=F, fix.mu=T, cond=cond, model = models[res_l]),
          silent = T
        ))
      } else{
        name <- paste0(models[res_l],"c")
        suppressWarnings(test <- try(
          treei_BCST <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                          backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                          f.lamb = function(x,y){y}, f.mu = function(x,y){0},
                                          lamb_par = list_param[[p]][1], mu_par = NULL,
                                          cst.lamb=T, cst.mu=T, expo.lamb=F, expo.mu=F, fix.mu=T, cond=cond, model = models[res_l],
                                          n.max = n.max, rate.max = rate.max),
          silent = T
        ))
      }
      p = p + 1
    }    
    
    if(verbose == T){cat("   ", name, " \t \t AICc =", treei_BCST$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda")] <- c(name, 1, c(treei_BCST$LH,treei_BCST$aicc, treei_BCST$lamb_par))
    res_l <- res_l + 1
    
  }
  if(treei_BCST$lamb_par > list_param[["improved"]][3]){
    list_param[["improved"]][1] <- treei_BCST$lamb_par
  }
  
  if("BCST_DCST" %in% models){
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        test <- tryCatch(
          treei_BCST_DCST <- fit_bd_backbone(phylo, tot_time, f=f,
                                             backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                             f.lamb = function(x,y){y}, f.mu = function(x,y){y},
                                             lamb_par = list_param[[p]][1], mu_par = list_param[[p]][3],
                                             cst.lamb=T, cst.mu=T, expo.lamb=F, expo.mu=F, fix.mu=F, cond=cond, model = models[res_l])
        )
      } else {
        name <- paste0(models[res_l],"c")
        test <- try(
          treei_BCST_DCST <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                               backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                               f.lamb = function(x,y){y}, f.mu = function(x,y){y},
                                               lamb_par = list_param[[p]][1], mu_par = list_param[[p]][3],
                                               cst.lamb=T, cst.mu=T, expo.lamb=F, expo.mu=F, fix.mu=F, cond=cond, model = models[res_l],
                                               n.max = n.max, rate.max = rate.max),
          silent = T
        )
      }
      p = p + 1
    }
    
    if(verbose == T){cat("   ", name, "\t \t AICc =", treei_BCST_DCST$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda", "Mu")] <- c(name, 2, c(treei_BCST_DCST$LH, treei_BCST_DCST$aicc, treei_BCST_DCST$lamb_par, treei_BCST_DCST$mu_par))
    res_l <- res_l + 1
    
  }
  
  if("BVAR" %in% models){
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        test <- try(
          treei_BVAR <- fit_bd_backbone(phylo, tot_time, f=f,
                                        backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                        f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){0},
                                        lamb_par = list_param[[p]][c(1,2)], mu_par = NULL,
                                        cst.lamb=F, cst.mu=T, expo.lamb=T, expo.mu=F, fix.mu=T, cond=cond, model = models[res_l]),
          silent = T
        )
        
        
      } else {
        name <- paste0(models[res_l],"c")
        test <- try(
          treei_BVAR <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                          backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                          f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){0},
                                          lamb_par = list_param[[p]][c(1,2)], mu_par = NULL,
                                          cst.lamb=F, cst.mu=T, expo.lamb=T, expo.mu=F, fix.mu=T, cond=cond, model = models[res_l],
                                          n.max = n.max, rate.max = rate.max),
          silent = T
        )
      }
      
      p = p + 1 
    }
    
    
    if(verbose == T){cat("   ",name," \t \t AICc =", treei_BVAR$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda", "Alpha")] <- c(name, 2, c(treei_BVAR$LH, treei_BVAR$aicc, treei_BVAR$lamb_par))
    res_l <- res_l + 1
  }
  
  list_param[["improved"]][c(1,2)] <- treei_BVAR$lamb_par
  
  if("BVAR_DCST" %in% models){
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        test <- try(
          treei_BVAR_DCST <- fit_bd_backbone(phylo, tot_time, f=f,
                                             backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                             f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){y},
                                             lamb_par = list_param[[p]][c(1,2)], mu_par = list_param[[p]][3],
                                             cst.lamb=F, cst.mu=T, expo.lamb=T, expo.mu=F,fix.mu=F,cond=cond, model = models[res_l]),
          silent = T
        )
      } else {
        name <- paste0(models[res_l],"c")
        test <- try(
          treei_BVAR_DCST <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                               backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                               f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){y},
                                               lamb_par = list_param[[p]][c(1,2)], mu_par = list_param[[p]][3],
                                               cst.lamb=F, cst.mu=T, expo.lamb=T, expo.mu=F,fix.mu=F,cond=cond, model = models[res_l],
                                               n.max = n.max, rate.max = rate.max),
          silent = T
        )
      }
      
      p = p + 1
    }
    
    if(verbose == T){cat("   ",name,"\t \t AICc =", treei_BVAR_DCST$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda", "Alpha", "Mu")] <- c(name, 3, c(treei_BVAR_DCST$LH, treei_BVAR_DCST$aicc, treei_BVAR_DCST$lamb_par, treei_BVAR_DCST$mu_par))
    res_l <- res_l + 1
  }
  
  if("BCST_DVAR" %in% models){
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        test <- try(
          treei_BCST_DVAR <- fit_bd_backbone(phylo, tot_time, f=f,
                                             backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                             f.lamb = function(x,y){y}, f.mu = function(x,y){y[1]*exp(y[2]*x)},
                                             lamb_par = list_param[[p]][c(1)], mu_par = list_param[[p]][c(3,4)],
                                             cst.lamb = T, cst.mu=F, expo.lamb=F, expo.mu=T, fix.mu=F, cond=cond, model = models[res_l]),
          silent = T
        )
      } else {
        name <- paste0(models[res_l],"c")
        test <- try(
          treei_BCST_DVAR <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                               backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                               f.lamb = function(x,y){y}, f.mu = function(x,y){y[1]*exp(y[2]*x)},
                                               lamb_par = list_param[[p]][c(1)], mu_par = list_param[[p]][c(3,4)],
                                               cst.lamb = T, cst.mu=F, expo.lamb=F, expo.mu=T, fix.mu=F, cond=cond, model = models[res_l],
                                               n.max = n.max, rate.max = rate.max),
          silent = T
        )
      }
      
      p = p + 1
    }
    
    if(verbose == T){cat("   ",name,"\t \t AICc =", treei_BCST_DVAR$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda", "Mu", "Beta")] <- c(name, 3, c(treei_BCST_DVAR$LH, treei_BCST_DVAR$aicc, treei_BCST_DVAR$lamb_par, treei_BCST_DVAR$mu_par))
    res_l <- res_l + 1
  }
  
  if("BVAR_DVAR" %in% models){
    p = 1
    test <- ""
    class(test) <- "try-error"
    while(class(test) == "try-error"){
      if(is.null(n.max) & is.null(rate.max)){
        name <- models[res_l]
        test <- try(
          treei_BVAR_DVAR <- fit_bd_backbone(phylo, tot_time, f=f,
                                             backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                             f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){y[1]*exp(y[2]*x)},
                                             lamb_par = list_param[[p]][c(1,2)], mu_par = list_param[[p]][c(3,4)],
                                             cst.lamb=F, cst.mu=F, expo.lamb=T, expo.mu=T, fix.mu=F, cond=cond, model = models[res_l]),
          silent = T
        )
      } else {
        name <- paste0(models[res_l],"c")
        test <- try(
          treei_BVAR_DVAR <- fit_bd_backbone_c(phylo, tot_time, f=f,
                                               backbone = backbone, spec_times = spec_times, branch_times = branch_times,
                                               f.lamb = function(x,y){y[1]*exp(y[2]*x)}, f.mu = function(x,y){y[1]*exp(y[2]*x)},
                                               lamb_par = list_param[[p]][c(1,2)], mu_par = list_param[[p]][c(3,4)],
                                               cst.lamb=F, cst.mu=F, expo.lamb=T, expo.mu=T, fix.mu=F, cond=cond, model = models[res_l],
                                               n.max = n.max, rate.max = rate.max),
          silent = T
        )
      }
      p = p + 1
    }
    
    if(verbose == T){cat("   ",name,"\t \t AICc =", treei_BVAR_DVAR$aicc,"\n")}
    results[res_l, colnames(results) %in% c("Models","Parameters","logL","AICc","Lambda", "Alpha", "Mu", "Beta")] <- c(name, 4, c(treei_BVAR_DVAR$LH, treei_BVAR_DVAR$aicc, treei_BVAR_DVAR$lamb_par, treei_BVAR_DVAR$mu_par))
    res_l <- res_l + 1
  }
  
  results[,-1] <- apply(results[,-1], 2, as.numeric)
  return(results)
}
