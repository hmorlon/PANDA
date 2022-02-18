remove.model <- function(shift.res, model){
  
  if(!model %in% c("BCST", "BCST_DCST", "BVAR", "BVAR_DCST", "BCST_DVAR", "BVAR_DVAR")){
    stop("Model name is incorrect. Should be the same name as in the function shift.estimates.")
  }
  
  whole <- shift.res$whole_tree
  whole <- whole[!whole$Models %in% c(model, paste0(model, "c")),]
  whole$delta_AICc <- whole$AICc - min(whole$AICc)
  whole <- whole[order(whole$AICc),]
  
  sub_list <- shift.res$subclades
  bck_list <- shift.res$backbones
  total2 <- as.data.frame(matrix(nrow = nrow(shift.res$total), ncol = ncol(shift.res$total)-1))
  names(total2) <- c("Combination", "Parameters", "logL", "AICc")
  
  sub_list <- lapply(seq_along(sub_list), function(i){
    sub_list[[i]] <- sub_list[[i]][!sub_list[[i]]$Models %in% c(model, paste0(model, "c")),]
    sub_list[[i]]$delta_AICc <- sub_list[[i]]$AICc - min(sub_list[[i]]$AICc)
    sub_list[[i]] <- sub_list[[i]][order(sub_list[[i]]$AICc), ]
  })
  
  names(sub_list) <- names(shift.res$subclades)
  
  # Multibackbone
  
  #is(sapply(bck_list, function(x) sapply(x, length)))[1] == "list"
  
  for(i in 1:length(bck_list)){
    for(j in 1:length(bck_list[[i]])){
      bck_list[[i]][[j]] <- bck_list[[i]][[j]][!bck_list[[i]][[j]]$Models %in% c(model, paste0(model, "c")),]
      bck_list[[i]][[j]]$delta_AICc <- bck_list[[i]][[j]]$AICc - min(bck_list[[i]][[j]]$AICc)
      bck_list[[i]][[j]] <- bck_list[[i]][[j]][order(bck_list[[i]][[j]]$AICc),]
    }
  }
  
  # summing up parts
  
  for(i in seq_along(bck_list)){
    comb <- unlist(strsplit(names(bck_list[i]), split = "/"))
    comb.sub <- comb[1]
    
    if(length(comb) == 2){
      comb.bck <- sapply(strsplit(names(bck_list[i]), split = "/"), "[[", 2) 
    }else{
      comb.bck <- ""
    }
    
    parts <- names(bck_list[[i]])
    best_parts <- lapply(bck_list[[i]], function(x) x[1, c("Parameters", "logL", "AICc")])
    
    # backbone parts
    total2[i, c("Parameters", "logL", "AICc")] <- apply(do.call(rbind, best_parts), 2, sum)
    total2$Combination[i] <- paste(comb.sub, comb.bck, sep = "/")
    
    # subclade parts
    comb.sub_split <- strsplit(comb.sub, split = "[.]")[[1]]
    for(n in 1:length(comb.sub_split)){
      total2[i, c("Parameters", "logL", "AICc")] <- total2[i, c("Parameters", "logL", "AICc")] +
        sub_list[[comb.sub_split[n]]][1, c("Parameters", "logL", "AICc")]
    }
  }
  
  total2$Combination[nrow(total2)] <- "whole_tree"
  total2[nrow(total2), c("Parameters", "logL", "AICc")] <- shift.res$whole_tree[shift.res$whole_tree$AICc == min(shift.res$whole_tree$AICc), c("Parameters", "logL", "AICc")]
  total2$delta_AICc <- total2$AICc - min(total2$AICc)
  total2 <- total2[order(total2$AICc),]
  row.names(total2) <- NULL
  
  shift.res2 <- list(whole_tree = whole, subclades = sub_list, backbones = bck_list, total = total2)
  
  return(shift.res2)
}
