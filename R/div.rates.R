div.rates <- function(phylo, shift.res, combi = 1, part = "backbone",
                      backbone.option = "crown.shift"){
  
  # Checking arguments ####
  # phylo
  if(!inherits(phylo, "phylo")){
    stop("object \"phylo\" is not of class \"phylo\"")
  }
  # shift.res
  if(!is(shift.res)[1] == "list" | any(names(shift.res) != c("whole_tree", "subclades", "backbones", "total"))){
    stop("object \"shift.res\" might be incorrect.")
  }
  if(!is.numeric(combi)){
    stop("object \"combi\" should be numeric.")
  }
  if(!part %in% c("backbone", "subclades", "all")){
    stop("object \"part\" might be incorrect.")
  }
  
  # Script ####
  
  # Subclades ####
  
  best_subclades_df <- do.call(rbind.data.frame, lapply(shift.res$subclades, function(x) x[1,]))
  best_subclades_df$Clades <- row.names(best_subclades_df)
  row.names(best_subclades_df) <- NULL
  best_subclades_df <- best_subclades_df[,c(10,1:8)]
  
  comb <- shift.res$total$Combination[combi]
  if(comb == "whole_tree"){
    
    tot_time <- max(branching.times(phylo))
    time.seq <- c(tot_time, seq(floor(tot_time),0,by=-1))
    
    rate_data <- shift.res$whole_tree[shift.res$whole_tree$AICc == min(shift.res$whole_tree$AICc),]
    model <- as.character(rate_data$Models)
    
    model <- rate_data$Models
    
    rate_df <- matrix(NA,2, length(time.seq))
    row.names(rate_df) <- c("Speciation", "Extinction")
    
    if (grepl("BCST", model)){
      rate_ext <- rep(NA, length(time.seq))
    }
    
    if (grepl("BCST_DCST", model)){
      rate_spec <- rep(rate_data$Lambda, length(time.seq))
      rate_ext <- rep(rate_data$Mu, length(time.seq))
    }
    
    if (grepl("BVAR", model)){
      rate_spec <- rate_data$Lambda *exp(rate_data$Alpha * time.seq)
      rate_ext <- rep(NA, length(time.seq))
    }
    
    if (grepl("BVAR_DCST", model)){
      rate_spec <- rate_data$Lambda *exp(rate_data$Alpha * time.seq)
      rate_ext <- rep(rate_data$Mu, length(time.seq))
    }
    
    if (grepl("BCST_DVAR", model)){
      rate_spec <- rep(rate_data$Lambda, length(time.seq))
      rate_ext <- rate_data$Mu *exp(rate_data$Beta * time.seq)
    }
    
    if (grepl("BVAR_DVAR", model)){
      rate_spec <- rate_data$Lambda *exp(rate_data$Alpha * time.seq)
      rate_ext <- rate_data$Mu *exp(rate_data$Beta * time.seq)
    }
    rate_df2 <- rbind(rate_spec, rate_ext)
    
    rate_df[,1:ncol(rate_df2)] <- rate_df2[,ncol(rate_df2):1]
    
    rownames(rate_df) <- c("Speciation", "Extinction")
    
    globalrate <- list(rate_df[,ncol(rate_df):1])
    names(globalrate) <- "whole_tree"
    
  } else {
    if(length(grep("/", comb)) == 1){
      if(length(strsplit(comb, "/")[[1]]) > 1){
        comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
        comb.bck <- strsplit(sapply(strsplit(comb, "/"), "[[", 2), "[.]")[[1]]
      } else{
        comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
        comb.bck <- NULL
      }
    } else {
      comb.sub <- strsplit(sapply(strsplit(comb, "/"), "[[", 1), "[.]")[[1]]
      comb.bck <- NULL
    }
    
    best_subclades_df_combi <- best_subclades_df[best_subclades_df$Clades %in% as.numeric(comb.sub),]
    names(best_subclades_df_combi)[1] <- "Parts"
    tot_time <- max(branching.times(phylo))
    
    time.seq <- c(tot_time, seq(floor(tot_time),0,by=-1))
    
    if(backbone.option == "stem.shift"){
      parental_nodes <- Ancestors(phylo, as.numeric(best_subclades_df_combi$Parts), type = "parent")
      tot_time2 <- as.list(branching.times(phylo)[as.character(parental_nodes)])
    } else {
      tot_time2 <- as.list(branching.times(phylo)[best_subclades_df_combi$Parts])
    }
    
    # Backbones #### 
    
    best_backbones <- shift.res$backbones[paste(paste(comb.sub, collapse = "."), paste(comb.bck, collapse = "."), sep = "/")][[1]]
    best_backbones_df <- do.call(rbind.data.frame, lapply(best_backbones, function(x) x[1,]))
    best_backbones_df$Parts <- row.names(best_backbones_df)
    row.names(best_backbones_df) <- NULL
    best_backbones_df <- best_backbones_df[,c(10,1:8)]
    
    
    if(part == "backbone"){
      
      rate_data <- best_backbones_df
      if(!is.null(comb.bck)){
        bck_names <- gsub("_sub", "", best_backbones_df$Parts)
        bck_names <- bck_names[-length(bck_names)]
        time_data <- c(branching.times(phylo)[bck_names], tot_time)
      } else {
        time_data <- tot_time
      }
      
      globalrate <- rep(list(NULL), nrow(rate_data))
      names(globalrate) <- c(best_backbones_df$Parts)
    }
    
    if(part == "subclades"){
      rate_data <- best_subclades_df_combi
      time_data <- unlist(tot_time2)
      globalrate <- rep(list(NULL), nrow(rate_data))
      names(globalrate) <- c(best_subclades_df_combi$Parts)
    }
    
    if(part == "all"){
      
      if(!is.null(comb.bck)){
        bck_names <- gsub("_sub", "", best_backbones_df$Parts)
        bck_names <- bck_names[-length(bck_names)]
        time_data <- c(unlist(tot_time2), branching.times(phylo)[bck_names], tot_time)
      } else {
        time_data <- c(unlist(tot_time2), tot_time)
      }
      
      rate_data <- rbind(best_subclades_df_combi, best_backbones_df)
      globalrate <- rep(list(NULL), nrow(rate_data))
      names(globalrate) <- c(best_subclades_df_combi$Parts, best_backbones_df$Parts)
    }
    
    for(r in 1:length(globalrate)){
      
      model <- as.character(rate_data$Models[r])
      agei <- time_data[r]
      time.seq_r <- c(tot_time, seq(floor(agei),0,by=-1))
      #time.seq_r <- unlist(ifelse(round(agei) == floor(agei), list(seq(floor(agei),0,by=-1)), list(c(agei, seq(floor(agei),0,by=-1)))))
      
      rate_df <- matrix(NA,2, length(time.seq))
      row.names(rate_df) <- c("Speciation", "Extinction")
      
      if (grepl("BCST", model)){
        rate_spec <- rep(rate_data$Lambda[r], length(time.seq_r))
        rate_ext <- rep(NA, length(time.seq_r))
      }
      
      if (grepl("BCST_DCST", model)){
        rate_spec <- rep(rate_data$Lambda[r], length(time.seq_r))
        rate_ext <- rep(rate_data$Mu[r], length(time.seq_r))
      }
      
      if (grepl("BVAR", model)){
        rate_spec <- rate_data$Lambda[r] *exp(rate_data$Alpha[r] * time.seq_r)
        rate_ext <- rep(NA, length(time.seq_r))
      }
      
      if (grepl("BVAR_DCST", model)){
        rate_spec <- rate_data$Lambda[r] *exp(rate_data$Alpha[r] * time.seq_r)
        rate_ext <- rep(rate_data$Mu[r], length(time.seq_r))
      }
      
      if (grepl("BCST_DVAR", model)){
        rate_spec <- rep(rate_data$Lambda[r], length(time.seq_r))
        rate_ext <- rate_data$Mu[r] *exp(rate_data$Beta[r] * time.seq_r)
      }
      
      if (grepl("BVAR_DVAR", model)){
        rate_spec <- rate_data$Lambda[r] *exp(rate_data$Alpha[r] * time.seq_r)
        rate_ext <- rate_data$Mu[r] *exp(rate_data$Beta[r] * time.seq_r)
      }
      rate_df2 <- rbind(rate_spec, rate_ext)
      
      rate_df[,1:ncol(rate_df2)] <- rate_df2[,ncol(rate_df2):1]
      
      rownames(rate_df) <- c("Speciation", "Extinction")
      
      globalrate[[r]] <- rate_df[,ncol(rate_df):1]
    }
  }
  
  return(globalrate)
}
