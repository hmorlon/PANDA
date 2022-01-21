extract.clade.ln<-function (phy, node, root.edge = 0) 
{
    Ntip <- length(phy$tip.label)
    ROOT <- Ntip + 1
    Nedge <- dim(phy$edge)[1]
    wbl <- !is.null(phy$edge.length)
    if (length(node) > 1) {
        node <- node[1]
        warning("only the first value of 'node' has been considered")
    }
    if (is.character(node)) {
        if (is.null(phy$node.label)) 
            stop("the tree has no node labels")
        node <- which(phy$node.label %in% node) + Ntip
    }
    if (node <= Ntip) 
        stop("node number must be greater than the number of tips")
    if (node == ROOT) 
        return(phy)
    	
    	
   	keep_nodes<-c(node)
   	
   	while(1)
   	{
   			keep_nodes1<-unique(c(keep_nodes,phy$edge[which(phy$edge[,1] %in% keep_nodes),2]))
   			if (length(keep_nodes1)==length(keep_nodes))
   			break
   			else
   			keep_nodes<-keep_nodes1
   			
   		}	
   	
    #print(keep_nodes)
    	
    keep<-which(phy$edge[,1] %in% keep_nodes)
   	
    phy$edge <- phy$edge[keep, ]
    if (wbl) 
        phy$edge.length <- phy$edge.length[keep]
    TIPS <- phy$edge[, 2] <= Ntip
    tip <- phy$edge[TIPS, 2]
    phy$tip.label <- phy$tip.label[tip]
    phy$edge[TIPS, 2] <- order(tip)
    if (!is.null(phy$node.label)) 
        phy$node.label <- phy$node.label[sort(unique(phy$edge[, 
            1])) - Ntip]
    Ntip <- length(phy$tip.label)
    phy$Nnode <- dim(phy$edge)[1] - Ntip + 1L
    newNb <- integer(Ntip + phy$Nnode)
    newNb[node] <- Ntip + 1L
    sndcol <- phy$edge[, 2] > Ntip
    phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]] <- (Ntip + 
        2):(Ntip + phy$Nnode)
    phy$edge[, 1] <- newNb[phy$edge[, 1]]
    phy
}
