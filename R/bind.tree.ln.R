bind.tree.ln<-function (x, y, where = "root", position = 0) 
# graft tree y on tree x at the node given by where, and distance from this node given by position. The distance is calculated from the node going backwards in time
{
	source("vect.from.string.r")
	source("string.from.vect.r")
	
    nb.tip <- length(x$tip.label)
    nb.node <- x$Nnode
    ROOT <- nb.tip + 1
    if (where == 0 || where == "root") 
        where <- ROOT
    if (position < 0) 
        position <- 0
    if (where > nb.tip + nb.node) 
        stop("node number out of range for tree 'x'")
    nb.edge <- dim(x$edge)[1]
    yHasNoRootEdge <- is.null(y$root.edge)
    xHasNoRootEdge <- is.null(x$root.edge)
    wbl <- TRUE
    noblx <- is.null(x$edge.length)
    nobly <- is.null(y$edge.length)
    if (noblx && nobly) 
        wbl <- FALSE
    if (xor(noblx, nobly)) {
        if (nobly) 
          	 x$edge.length <- NULL
        else y$edge.length <- NULL
        wbl <- FALSE
        warning("one tree has no branch lengths, they will be ignored")
    }
    if (where <= nb.tip) {
        Tip.Label.where <- x$tip.label[where]
        x$tip.label[where] <- "TheTipWhereToGraftY"
    }
    if (where > ROOT) {
        xHasNoNodeLabel <- TRUE
        if (is.null(x$node.label)) {
            x$node.label <- paste("NODE", 1:nb.node, sep = "")
            x$node.label[where - nb.tip] <- "TheNodeWhereToGraftY"
        }
        else {
            Node.Label.where <- x$node.label[where - nb.tip]
            x$node.label[where - nb.tip] <- "TheNodeWhereToGraftY"
            xHasNoNodeLabel <- FALSE
        }
    }
    if (wbl) {
        if (where == ROOT) {
            new.root <- max(x$root.edge - position,0)
            x$root.edge <- position
        }
        else {
            i <- which(x$edge[, 2] == where)
            if (x$edge.length[i] < position) 
                stop("argument 'position' is larger than the specified edge.")
            x$edge.length[i] <- x$edge.length[i] - position
        }
        if (yHasNoRootEdge) 
            y$root.edge <- 0
    }
    X <- write.tree(x)
    #print(X)
    Y <- write.tree(y)
    Y <- substr(Y, 1, nchar(Y) - 1)
    if (where <= nb.tip) {
        if (position) 
            {X <- gsub("TheTipWhereToGraftY", paste("(", "TheTipWhereToGraftY", ":", position,
                ",", Y, ")", sep = ""), X)
                #print(X)
                }
        else 
        stop("the code is not well implemented for this case")
        #X <- gsub("TheTipWhereToGraftY", Y, X)
    }
    else if (where == ROOT) {
    	X <- substr(X, 1, nchar(X) - 1)       
        X <- paste("(", X, ",", Y, ")", ":", new.root, ";", sep = "")
    }
    else if (where > ROOT) { 
    		#print(X)	
        	X.vect<-vect.from.string(X)
        	X.length<-length(X.vect)
        	#print(X.vect)
        	TheNodeWhereToGraftY.position<-which(X.vect=="TheNodeWhereToGraftY")
        	character<-TheNodeWhereToGraftY.position-2
        	#print(X.vect[character])
        	#print(character)
        	count<-1
        	nbchar<-0
        	while(count>0)
        	{
        		nbchar<-(nbchar+1)
        		if (X.vect[character]==")")
        		{
        			count<-count+1
        			character<-character-1}
        		else if (X.vect[character]=="(")
        		{
        			count<-count-1
        			character<-character-1}
        		
        		else {character<-character-1}
        			}
        			#print(string.from.vect(X.vect[1:(character+nbchar)]))
        			#print(string.from.vect(X.vect[(TheNodeWhereToGraftY.position+3):X.length]))
        			X<-string.from.vect(c(X.vect[1:character],"(",X.vect[(character+1):(character+nbchar+1)],":",position,",",Y,")",":",x$edge.length[i],X.vect[(TheNodeWhereToGraftY.position+3):X.length],sep=""))
        			#print(X)

    }
    #print(count)
    phy <- read.tree(text = X)
    if (where <= nb.tip) 
        phy$tip.label[which(phy$tip.label == "TheTipWhereToGraftY")] <- Tip.Label.where
    if (where > ROOT) {
        if (xHasNoNodeLabel) 
            phy$node.label <- NULL
        else phy$node.label[which(phy$node.label == "TheNodeWhereToGraftY")] <- Node.Label.where
    }
    phy
}