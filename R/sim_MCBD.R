
sim_MCBD <- function (pars, root.value=0, age.max=50, step.size=0.01, bounds=c(-Inf,Inf), plot=T, ylims=NULL, full.sim=F){   

  lambda1 = pars[1]
  tau0 = pars[2]
  beta = pars[3]
  mu0 = pars[4]
  mubg = pars[5]
  mui0 = pars[6]
  muibg = pars[7]
  alpha1 = pars[8]
  alpha2 = pars[9]
  sig2 = pars[10]
  m = pars[11]
  
  s = sqrt(sig2 * step.size)
  m = m * step.size
  
  if (length(bounds)!=2){stop("Lower and upper bounds should be provided")}
  if (root.value < bounds[1] || root.value > bounds [2]) {print("Warning: root value outside boundaries; continuing simulation")}
  
  process_dead = F
  traits <- list(c(1,1,2,root.value),c(2,-1,1,root.value)) #lineage number, status (incipient=-1/good=1/extinct=-2), sister lineage, trait values 
  lineages <- rbind(c(1,0,0,-1,1,0,0),c(1,0,0,-1,-1,0,NA)) #parental node, descendant node, starting t, ending t (still alive=-1), status, sp completion/extinction t, sp completion t
  active_lineages = c(1,2)
  n_lineages = length(active_lineages)
  n_good = 2L
  step_count = 1L
  t = 0 + step.size
  
  #step by step simulation
  while ((age.max-t) > -step.size/2){#t must go one step beyond age.max to ensure simulation of last step. /2 is to avoid numerical precission issues
    mat_diff <- matrix(0, n_lineages, n_lineages, dimnames = list(active_lineages, active_lineages))
    
    #trait differences among lineages
    for (i in 1:n_lineages){
      mat_diff[i,] <- sapply(traits[active_lineages],function(x)(traits[[active_lineages[i]]][3+step_count] - x[3+step_count]))
    }
    
    diag(mat_diff) <- NA
    diff_sign <- sign(mat_diff)
    #if there are any equal lineages, assigns random and opposing repulsion directions
    if (any(diff_sign==0, na.rm = T)){
      diff_sign[lower.tri(diff_sign)] <- NA
      eq_ind <- which(diff_sign == 0, arr.ind = T, useNames = F)
      for (i in 1:nrow(eq_ind)){diff_sign[eq_ind[i,1],eq_ind[i,2]] <- sign(rnorm(1))}
      diff_sign[lower.tri(diff_sign)] <- -diff_sign[upper.tri(diff_sign)]
      }
    
    #traits loop
    diff_me = list()
    for (i in 1:length(active_lineages)){
      signs_me <- diff_sign[i,-i]
      diff_me[[active_lineages[i]]] <- mat_diff[i,-i] #saves diff vector to use later for extinction rates
      dist_bounds <- traits[[active_lineages[i]]][3+step_count] - bounds #distance to bounds
      bound_effect <- 3 * sum(sign(dist_bounds)*exp(-2*(dist_bounds)^2)) #repulsion of bounds
      #update trait value
      traits[[active_lineages[i]]][4+step_count] <- traits[[active_lineages[i]]][3+step_count] + 
        alpha2 * m * sum(signs_me * exp(-alpha2 * (diff_me[[active_lineages[i]]])^2)) +
        rnorm(1, 0, s) + bound_effect 
    }
    
    dead_lin = NULL
    born_lin = NULL
    #lineages branching loop
    for (i in active_lineages){
      #if lineage is incipient
      if (lineages[i,5] == -1){
        #captures the case where sister lineage speciates again before finishing another i_sp  
        if (is.na(traits[[traits[[i]][3]]][[3+step_count]])){
          diff_trait_isp = abs(traits[[i]][3+step_count] - traits[[traits[[i]][3]]][[4+step_count]])
        }
        else{diff_trait_isp = abs(traits[[i]][3+step_count] - traits[[traits[[i]][3]]][[3+step_count]])}
        #sp completion rate depends on distance w/parent
        lambda2 = tau0 * exp(beta * (diff_trait_isp)^2)
        #extinction rate depends on distance with all lineages, same as trait evolution
        mui = alpha1 * mui0 * sum(exp(-alpha1 * (diff_me[[i]])^2)) + muibg 
        
        probs_i = c(lambda2/(mui+lambda2), mui/(mui+lambda2))
        if (runif(1) <= (lambda2 + mui) * step.size) {
          event <- sample (1:2, size = 1, prob = probs_i)
          #speciation completion
          if (event == 1){
            lineages[i,5:7] <- c(1, t, t)
            traits[[i]][2] <- 1
          }
          #extinction
          if (event == 2){
            lineages[i,c(4,5,6)] <- c(t,-2,t)
            traits[[i]][2] <- -2
            dead_lin <- c(dead_lin, i)
            # i sister has now new sister: it's most closely related living good sp (if any, otherwise takes closest incipient)
            i_sister <- traits[[i]][3]
            #to find it... take first ancestral node, any descendent living lineage? if fails, take second ancestral node, an so on. if many lineages, pick random.
            anc_node <- lineages[which(lineages[,2]==lineages[i_sister,1]),1]
            while (T){
              #collect all living tips from focal anc_node
              possible <- active_lineages[lineages[active_lineages,1]>=anc_node] #lineages that could be descendent of anc_node (i.e. younger)
              possible <- possible[possible!=i_sister] #not myself
              possible <- possible[lineages[possible,5]!=-2] #not extinct at t
              if(length(possible)==0){break}
              #discard lineages that don't descend from anc_node
              for (lin in possible){
                node_parent <- lineages[lin,1]
                while(T){
                  if(node_parent==anc_node){break}
                  else if(node_parent<anc_node){possible<-possible[possible!=lin]; break}
                  node_parent <- lineages[which(lineages[,2]==node_parent)]
                }
              }
              if (length(possible)>0){
                goods <- possible[lineages[possible,5]==1] 
                #if there are good lineages, choose goods first (randomly if there are many)
                if (length(goods)>0){traits[[i_sister]][3]<-goods[sample(1:length(goods),size = 1)]}
                else{traits[[i_sister]][3]<-possible[sample(1:length(possible),size = 1)]}
                break
              }#found!
              #otherwise move one node backwards and repeat
              anc_node <- lineages[which(lineages[,2]==anc_node),1]
            }
          }
        }
      }
      #if lineage is good
      if (lineages[i,5] == 1) {
        #extinction rate depends on distance with all lineages, same function as trait evolution
        mu = alpha1 * mu0 * sum(exp(-alpha1 * (diff_me[[i]])^2)) + mubg
        #captures the case when there's only one lineage alive & extinction is activated
        if (mu0!=0 & mu==0){
          mu = 0.02 * mu0 #when alone, a lineage has basal extinction rate (equal to having infinite distance with neighbors)
        }
        probs <- c(lambda1/(lambda1+mu), mu/(lambda1+mu))
        if (runif(1) <= (mu + lambda1) * step.size) {
          event <- sample(1:2, size = 1, prob = probs)
          #speciation
          if (event == 1){
            lineages[i,2] <- max(lineages[,1]) + 1
            lineages[i,4] <- t
            new_lin1 = c(lineages[i,2],0,t,-1,1, t,t)
            new_lin2 = c(lineages[i,2],0,t,-1,-1, NA, NA)
            lineages <- rbind(lineages, new_lin1, new_lin2)
            dead_lin <- c(dead_lin,i)
            born_lin <- c(born_lin, dim(lineages)[1]-1, dim(lineages)[1])
            traits[[dim(lineages)[1]-1]] <- c(dim(lineages)[1]-1,1,dim(lineages)[1],rep(NA, times = step_count),traits[[i]][4+step_count])
            traits[[dim(lineages)[1]]] <- c(dim(lineages)[1],-1,dim(lineages)[1]-1,rep(NA, times = step_count),traits[[i]][4+step_count])
            #has already descendent incipient lineages?
            sisters <- which(unlist(lapply(traits,function(x)(x[3])))==i) #get lineages who have i as sister
            sisters <- sisters[lineages[sisters,2]==0] #keep only living 
            if (sum(sisters)>0){
              #update sister number of all i_sp lineages descending from this parental good lineage (or it's ancestors) with new id 
              for (j in sisters){
               traits[[j]][3] <- nrow(lineages)-1
              }
            }
          }
          #extinction
          if (event == 2){
            lineages[i,c(4,5,6)] <- c(t,-2,t)
            traits[[i]][2] <- -2
            dead_lin <- c(dead_lin, i)
            #if sister is incipient & alive, becomes good 
            i_sister <- traits[[i]][3]
            if (lineages[i_sister,5] == -1){lineages[i_sister,c(5,6,7)] <- c(1,t,t); traits[[i_sister]][2]<-1}
            #i_sister inherits older sisters from the extinct
            sisters <- which(unlist(lapply(traits,function(x)(x[3])))==i) #get older lineages who have i as sister
            #discard current sister, extinct, and ancestral lineages
            sisters <- sisters[sisters!=i_sister] 
            sisters <- sisters[which(lineages[sisters,5]!=-2)]
            sisters <- sisters[which(lineages[sisters,2]==0)]
            for (j in sisters){traits[[j]][3] <- i_sister}
            #i_sister has new sister, which is the youngest of the older non-extinct branches of i (if any left)
            if (sum(sisters)>0){traits[[i_sister]][3] <- max(sisters)}
            }
          }
      }
      }
    
    if (!is.null(dead_lin)){active_lineages <- c(active_lineages[!(active_lineages%in%dead_lin)], born_lin)}
    
    step_count = step_count + 1L
    t = t + step.size
    cat("\r","time:", t)
    n_lineages = length(active_lineages)
    
    if (n_lineages == 0) {print ("process died")
        process_dead = T
        break()
      }
    }
  t = t - step.size
  
  
  ##BUILD TREES & GET TIPS TRAIT VECTORS
  #complete process tree
  row.names(lineages) <- NULL
  colnames(lineages) <- NULL
  edges_mat <- lineages[,1:2]
  active_lineages <- sort(c(active_lineages, which(lineages[,5]==-2)))
  n_tips = length(active_lineages)
  edges_mat[,1] <- edges_mat[,1] + n_tips
  edges_mat[,2][which(edges_mat[,2]!=0)] <- edges_mat[,2][which(edges_mat[,2]!=0)] + n_tips
  edges_mat[,2][which(edges_mat[,2]==0)] <- 1:n_tips
  edg1 <- as.integer(edges_mat[,1]); edg2 <- as.integer(edges_mat[,2]); edges_mat <- cbind(edg1,edg2) #to comply with phylo class
  dimnames(edges_mat) <- NULL
  lineages[,4][which(lineages[,4]==-1)] <- t
  lin_length <- round(lineages[,4] - lineages[,3],2)
  tree <- list(edge = edges_mat, edge.length = lin_length, Nnode = (n_tips-1), tip.label = paste("t", as.character(1:n_tips), sep="") )
  class(tree) <- "phylo"
  #traits
  tip_values = NULL
  tip_values <- sapply(traits[active_lineages], function(x)(x[length(x)]))
  names(tip_values) <- tree$tip.label
  ###
  #good species tree (extant & extinct)
  isp_todrop <- c(which(lineages[,5]==-1),which(lineages[,5]==-2&is.na(lineages[,7]))) 
  tip_ids <- edges_mat[isp_todrop,2]
  tree_gsp_fossil <- drop.tip(tree, tip = tip_ids)
  #traits
  tip_values_gsp_fossil <- tip_values[names(tip_values) %in% tree_gsp_fossil$tip.label]
  
  if(process_dead == T){
    tree_gsp_extant <- "process died"
    tip_values_gsp_extant <- "process died"
  }
  else{
    #reconstructed good species tree
    extinct_todrop <- c(which(lineages[,5]==-1),which(lineages[,5]==-2))
    tip_ext_ids <- edges_mat[extinct_todrop,2]
    tree_gsp_extant <- drop.tip(tree, tip = tip_ext_ids)
    #traits
    tip_values_gsp_extant <- tip_values[names(tip_values) %in% tree_gsp_extant$tip.label]
  }
  
  #to get internal order of trees as commonly expected
  tree <- read.tree(text=write.tree(tree))
  tree_gsp_fossil <- read.tree(text=write.tree(tree_gsp_fossil))
  tree_gsp_extant <- read.tree(text=write.tree(tree_gsp_extant))
  
  
  res <- list(all=list(tree=tree, tips=tip_values),
              gsp_fossil=list(tree=tree_gsp_fossil,tips=tip_values_gsp_fossil), 
              gsp_extant=list(tree=tree_gsp_extant,tips=tip_values_gsp_extant))
  
  if (full.sim == T){
    res$all$trait_mat <- traits
    res$all$lin_mat <- lineages
  }
  
  if (plot){
    plotSimu <- function(traitmat, linmat, step_size, ylims=NULL){
      #plots a simulation, with incipient lineages in red, good in black
      steps <- max(linmat[,4])/step_size
      max_trait <- NULL
      min_trait <- NULL
      for (i in 1:nrow(linmat)){#find extremes for plot limits
        max_trait <- c(max_trait, max(traitmat[[i]][-(1:3)], na.rm = T))
        min_trait <- c(min_trait, min(traitmat[[i]][-(1:3)], na.rm = T))
      }
      if (is.null(ylims)){
        plot(1, type="n",xlim=c(1,steps), ylim = c(min(min_trait),max(max_trait)))
      }
      else{plot(1, type="n",xlim=c(1,steps), ylim = ylims)}
      completion <- linmat[,6]
      #handle and plot each possibility
      for (i in 1:length(completion)){
        if (is.na(completion[i])){
          #extinct incipient lineages
          lines(x = (4:length(traitmat[[i]])), traitmat[[i]][4:length(traitmat[[i]])], col="red", cex=0.5)
        } 
        else if (linmat[i,5] == -2 && !is.na(linmat[i,7]) && completion[i]!= linmat[i,7]){
          #extinct good lineages
          lines(x = 4:((linmat[i,7]/step_size)+4), y = traitmat[[i]][4:((linmat[i,7]/step_size)+4)], col="red", cex=0.5)
          lines(x=((linmat[i,7]/step_size)+4):length(traitmat[[i]]), y=traitmat[[i]][((linmat[i,7]/step_size)+4):length(traitmat[[i]])], col="black", cex=0.5)
        }
        else{
          #living good and incipient lineages
          lines(x = 4:((completion[i]/step_size)+4), y = traitmat[[i]][4:((completion[i]/step_size)+4)], col="red", cex=0.5)
          if(length(traitmat[[i]])>(round(completion[i]/step_size)+4)){
            lines(x=((completion[i]/step_size)+4):length(traitmat[[i]]), y=traitmat[[i]][((completion[i]/step_size)+4):length(traitmat[[i]])], col="black", cex=0.5)
          }
        }
      }
    }
    if(length(ylims)<2){ylims=NULL}
    plotSimu(traits,lineages, step.size, ylims)
  }
  
  return (res)
}
