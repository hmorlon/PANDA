##this version allows for a tree to not have sequentially numbered nodes (e.g., trees loaded with read.tree function)

require(geiger)
require(phytools)
require(deSolve)

.VCV.rescale.geog<-function(phylo,sigma,alpha,sterm,geography.object){
##these are the three parameters that the equations use to produce the variance-covariance matrix, and will ultimately be estimated with ML
parameters<-c(a=alpha,b=sigma,s=sterm) 
if(any(grepl("___",phylo$tip.label))){stop("VCV.rescale cannot handle trees with '___' in the tip labels")}

paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS


nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
for (i in 1:dim(phylo$edge)[1]){
	nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
}
nodeDist<-c(nodeDist,max(heights))
nodeDiff<-diff(nodeDist)
##label the branches for each segment of tree to be integrated and identify the node at which the branch terminates

if(sum(nodeDiff<0)>0){  ##this loop renumbers the nodes if trees nodes are not placed in sequential order
	node.order<-match(rank(heights[,1],ties.method="min"),seq(1, by = 2, len = phylo$Nnode))
	node.order<-node.order+length(phylo$tip.label)
	old.edge<-phylo$edge
	phylo$edge[,1]<-node.order
	for(j in 1:length(phylo$edge[,2])){
		if(phylo$edge[j,2]>length(phylo$tip.label)){
			#match number order in old edge
			#lookup value in new edge
			#replace with value
			phylo$edge[j,2]<-phylo$edge[,1][match(phylo$edge[j,2],old.edge[,1])]
			}
		}
	nodeDist<-vector()
	for (i in 1:dim(phylo$edge)[1]){
		nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
		}
	nodeDist<-c(nodeDist,max(heights))
	nodeDiff<-diff(nodeDist)
}
if(any(nodeDiff==0)){stop("VCV.rescale cannot handle trees with two or more nodes occurring at exactly the same time")}
if(any(nodeDiff<9e-6)){stop("The current tree has very small branch lengths (< 9e-6), consider multiplying all branch lengths to increase reliability of numerical integration")}
if(length(geography.matrix)!=phylo$Nnode){stop("The number of sympatry/allopatry matrices does not equal the number of time periods")}

mat<-matrix(nrow=0, ncol=3)
counter_three_letters <- 0
for(i in 1:phylo$Nnode){
	other<-phylo$edge[phylo$edge[,1]==i+length(phylo$tip.label), 2]
	for(b in other){
		int<-matrix(ncol=3)
		int[1]<-i+length(phylo$tip.label)
		if(b>length(phylo$tip.label)){
			counter_three_letters <- counter_three_letters + 1
			int[2]<-paste(".",THREELETTERS[counter_three_letters],sep="")
			int[3]<-b
			} else {
			int[2]<-phylo$tip.label[[b]]
			int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
			}
		mat<-rbind(mat,int)
		}
	}		
nat<-list()
for(i in 1:length(nodeDiff)){
	if(i==1){
	nat[[i]]<-list(mat[mat[,1]==(length(phylo$tip.label)+i),2])} else {
	IN<-vector()
	P<-mat[as.numeric(mat[,1])<=(length(phylo$tip.label)+i),c(2,3)]
	IN<-c(IN, P[P[,2]=="0",1],P[as.numeric(P[,2])>(length(phylo$tip.label)+i),1])
	nat[[i]]<-list(IN)
	}
	}	
for(i in 2:length(nodeDiff)){			##THIS LOOP checks for an error
	if(length(unlist(nat[[i]]))!=(length(unlist(nat[[i-1]]))+1)){
		print(paste("ERROR at node",i+length(phylo$tip.label)))
		}	
	}


#### NOW DEFINE ODEs for each interval and numerically integrate from root to tip, one interval at at a time ###
output<-list() ##initialize storage of results from ODEs 
env.results<-new.env(size = 1, parent = emptyenv())

for(i in 1:phylo$Nnode){ ##for each interval between branches...

##var.list is the list of terms for which there is a variance term at each step in time
##cov.list is the list of covariance terms
var.list<-unlist(nat[[i]])
cov.list<-apply(t(combn(var.list,2)),1,paste,collapse="___")
Cmat<-matrix(nrow=length(var.list),ncol=length(var.list)) ##"Cmat" is a lower triangle matrix with each of the variance and covariance terms, used to define equations below
diag(Cmat)<-var.list
Cmat[lower.tri(Cmat)]<-cov.list

##define ODEs for deSOLVE
	len <- length(var.list)
	dim <- len * (len + 1) / 2
	ind <- matrix(0, nrow = len, ncol=len, byrow=TRUE)
	indI <- vector(mode = "integer", length = dim)
	indJ <- vector(mode = "integer", length = dim)
	for(k in 1:len) {
	    ind[[k,k]] <- k
	    indI[[k]] <- k
	    indJ[[k]] <- k
	}
	counter <- len
	for(l in 1:(len-1)) {
	    for(k in (l+1):len) {
	        counter <- counter + 1
	        ind[[k,l]] <- counter
	        ind[[l,k]] <- counter
	        indI[[counter]] <- k
	        indJ[[counter]] <- l
	    }
	}
	coefAlpha <- (2*(len-1)/len * parameters['s'] + 2 * parameters['a'])
	coefBeta <- parameters['s'] / len
	rhs <- vector(mode = "numeric", length = dim)
	for(m in 1:len) {
	    rhs[m] <- parameters['b']
	}
	A <- diag(x = -coefAlpha, nrow = dim, ncol=dim)
	if(parameters['s']!=0){
	updated.alpha<-vector()	
	for(term in 1:length(indJ)){
		j<-indJ[term]
		int<-(((sum(geography.matrix[[i]][,j])-1)/sum(geography.matrix[[i]][,j]))/((len-1)/len)*(coefAlpha-2* parameters['a']))+2* parameters['a']
		updated.alpha<-c(updated.alpha, int)
		}	
	diag(A)<--updated.alpha
	}
	coefBetaV <- rep(coefBeta, times = len)
	for(term in 1:len){
		int<-coefBeta*(len/sum(geography.matrix[[i]][,term]))
		coefBetaV[term]<-int
		}	
	for(ind_jk in 1:dim) {
	    j <- indI[[ind_jk]]
	    k <- indJ[[ind_jk]]
	    coefBetaV[[k]] <- 0
	    #term<-ifelse(ind_jk%%len==0,len,ind_jk%%len)
	    if(geography.matrix[[i]][match(ind_jk,ind)]==0){A[ind_jk, ind[j, ] ] <- 0} else{ 
	    	##first make a list of match terms
	    	geog<-geography.matrix[[i]][match(ind[j,],ind)]
	    	A[ind_jk, ind[j, ] ] <- (A[ind_jk, ind[j, ] ] + coefBetaV)*geog
	    	} #should diagonal be preserved? does it matter?
	    coefBetaV[[k]] <- coefBeta*(len/sum(geography.matrix[[i]][,k]))
	    coefBetaV[[j]] <- 0
	    if(geography.matrix[[i]][match(ind_jk,ind)]==0){A[ind_jk, ind[k, ] ] <- 0} else{ 
	    	geog<-geography.matrix[[i]][match(ind[k,],ind)]
	    	A[ind_jk, ind[k, ] ] <- (A[ind_jk, ind[k, ] ] + coefBetaV)*geog
	    	} #should diagonal be preserved? does it matter?
 
	    coefBetaV[[j]] <- coefBeta*(len/sum(geography.matrix[[i]][,j]))
	}
	#diag(A)<-updated.alpha #whether this is included does affect result, but not sure how yet, should in theory contain the OU part of the model
	diag(A)<-sapply(diag(A),function(x) {if(x==0){x<--(2* parameters['a'])}else{x<-x}})
	
 # return the rate of change
 sttrt=paste("d",c(var.list,cov.list),sep="",collapse=",")
 sttat=paste("list(c(",sttrt,"))",sep="")

######	STARTING VALUES  ######
#This next block of code provides the starting values for the ODEs by identifying the proper value from the previous integration
if(i==1){ #if it is the first iteration, the starting values for all terms are 0
start=c(paste(c(var.list,cov.list),"=0",sep="",collapse=","))
state.new<-paste("c(",start,")",sep="")
state<-eval(parse(text=state.new))
} else {
#initialize starting value named vector
start=c(paste(c(var.list,cov.list),"=0",sep="",collapse=","))
state.new<-paste("c(",start,")",sep="")
state<-eval(parse(text=state.new))
termlist<-c(var.list,cov.list)
for(l in 1:length(termlist)){
term<-termlist[l]
if(exists(term, envir=env.results)){ #if the term to be numerically integrated was present in the previous generation(branch) of the tree
       state[l]<-get(term, envir=env.results)
	}  
else { #if term is not present in previous generation/branch of the tree, lookup which branch it descends from
	if(l <= length(var.list)){ ##this loop looks up values for variance terms using the "mat" matrix that is produced earlier
		prev.branch<-mat[mat[,3]==mat[mat[,2]==term,1],2]
		state[l]<-get(prev.branch, envir=env.results)
	}
	else{ ##for covariance terms that weren't present in previous generation/branch, it is necessary to decompose covariance term, look up terms appropriately, and re-assemble them in the correct order so they match the lower triangle format used in definition of terms
		firstterm<-strsplit(term,"___")[[1]][1]
		scondterm<-strsplit(term,"___")[[1]][2]
		if(exists(firstterm, envir=env.results) || exists(scondterm, envir=env.results)){
			if(exists(firstterm, envir=env.results)){
				prev.branch<-mat[mat[,3]==mat[mat[,2]==scondterm,1],2]
				prev.branch2<-firstterm} 
			else{  ##NOTE prev.branch2 just identifies which term was present previously
				prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
				prev.branch2<-scondterm
				}
		prev.term<-paste(prev.branch2,prev.branch,sep="___")
		if(! exists(prev.term, envir=env.results)){
			prev.term<-paste(prev.branch,prev.branch2,sep="___")
			}
		state[l]<-get(prev.term, envir=env.results)
		} 				
	else {
		##neither first nor second term was present in the previous generation/branch, so the starting value is a variance term from t-1
		prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
		state[l]<-get(prev.branch, envir=env.results)
	}
	}
	}
	}
	}

ou <- function(t, state, parameters) {
  dX <- A %*% state + rhs
  return (list(c(dX)))
}

##NOW, run numerical integration for given time step and append to a list
output<-ode(y=state,times=c(0,nodeDiff[i]),func=ou,parms=NULL)
colN <- colnames(output)
env.results<-new.env(size = length(colN), parent = emptyenv())
for (k in 2:length(colN)){
  assign(colN[[k]], output[[2,k]], envir = env.results)
}
}

Vou<-matrix(nrow=length(phylo$tip.label),ncol=length(phylo$tip.label))
rownames(Vou)<-unlist(nat[[phylo$Nnode]])
colnames(Vou)<-unlist(nat[[phylo$Nnode]])
diag(Vou)<-output[2,][2:(length(phylo$tip.label)+1)]
for(j in (length(phylo$tip.label)+2):length(output[2,])){

	string<-strsplit(names(output[2,j]),"___")
	string<-unlist(string)
	one<-match(string[1],colnames(Vou))
	two<-match(string[2],colnames(Vou))
	Vou[one,two]<-output[2,j]
	Vou[two,one]<-output[2,j]
}

return(Vou)
}

