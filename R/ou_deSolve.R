require(Matrix)
require(geiger)
require(phytools)
require(deSolve)
sigma=0.05
alpha=0.05
sterm=0.05

VCV.rescale<-function(phylo,sigma,alpha,sterm){
##these are the three parameters that the equations use to produce the variance-covariance matrix, and will ultimately be estimated with ML
parameters<-c(a=alpha,b=sigma,s=sterm) 

###FIRST make a vector with all unique three letter strings for naming internal branches, works up to 17576 tips (could be done once in session and then removed from function)

paste(rep(LETTERS,each=26),LETTERS,sep="")->TWOLETTERS
paste(rep(TWOLETTERS,each=26),LETTERS,sep="")->THREELETTERS


##define a few useful variables
##nodeDist, the distance between each inner node and the root
##nodeDiff, the distance between each inner node and the previous node (times for integration)

nodeDist<-vector(mode = "numeric", length = phylo$Nnode)
root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
for (i in 1:dim(phylo$edge)[1]){
	nodeDist[[phylo$edge[i, 1] - length(phylo$tip.label)]] <- heights[i]
}
nodeDist<-c(nodeDist,max(heights))
nodeDiff<-diff(nodeDist)
##label the branches for each segment of tree to be integrated and identify the node at which the branch terminates

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


##now come up with a list of branches that exist for each time segment (+1 at each branching in this version, which means tree can't have any polytomies or speciation events at exactly the same time)
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

require(Rcpp)
sourceCpp( "../src/fillIndices.cpp" )


for(i in 1:phylo$Nnode){ ##for each interval between branches...

##var.list is the list of terms for which there is a variance term at each step in time
##cov.list is the list of covariance terms
var.list<-unlist(nat[[i]])
cov.list<-apply(t(combn(var.list,2)),1,paste,collapse="_")

##define ODEs for deSOLVE
	##this loop defines the set of ODEs for the variance terms (stored in "var.list"), one for each "extant branch" of the phylo
	##note: if you want to visualize the output of this, replace the line beginning with 'eval(' with:
	##then call each term in var.list (i.e., in var.list) to see what the ODE for that term is
	len <- length(var.list)
	dim <- len * (len + 1) / 2
	coefAlpha <- (2*(len-1)/len * parameters['s'] + 2 * parameters['a'])
	coefBeta <- parameters['s'] / len
	rhs <- vector(mode = "numeric", length = dim)
	for(m in 1:len) {
	    rhs[m] <- parameters['b']
	}
	csrX <- numeric(len*len*(len-1))
	csrJ <- integer(len*len*(len-1))
	csrP <- integer(dim + 1)
	fillIndices(csrJ, csrP, csrX, len)

	spA <- sparseMatrix(j=csrJ, p=csrP, x=csrX)
	rm(csrJ, csrP, csrX)

	 
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
for(l in seq_along(termlist)){
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
		firstterm<-strsplit(term,"_")[[1]][1]
		scondterm<-strsplit(term,"_")[[1]][2]
		if(exists(firstterm, envir=env.results) || exists(scondterm, envir=env.results)){
			if(exists(firstterm, envir=env.results)){
				prev.branch<-mat[mat[,3]==mat[mat[,2]==scondterm,1],2]
				prev.branch2<-firstterm} 
			else{  ##NOTE prev.branch2 just identifies which term was present previously
				prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
				prev.branch2<-scondterm
				}
		prev.term<-paste(prev.branch2,prev.branch,sep="_")
		if(! exists(prev.term, envir=env.results)){
			prev.term<-paste(prev.branch,prev.branch2,sep="_")
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
  dX <- coefBeta * ((spA %*% state)@x) - coefAlpha * state + rhs
  return (list(c(dX)))
}

##NOW, run numerical integration for given time step and append to a list
output<-ode(y=state,times=c(0,nodeDiff[i]),func=ou,parms=NULL, method="adams")
colN <- colnames(output)
env.results<-new.env(size = length(colN), parent = emptyenv())
for (k in 2:length(colN)){
  assign(colN[[k]], output[[2,k]], envir = env.results)
}
}


#########	RETURN VCV MATRIX  	#########

Vou<-mrca(phylo)
diag(Vou)<-output[2,][2:(length(phylo$tip.label)+1)]
for(j in (length(phylo$tip.label)+2):length(output[2,])){

	string<-strsplit(names(output[2,j]),"_")
	string<-unlist(string)
	one<-match(string[1],colnames(Vou))
	two<-match(string[2],colnames(Vou))
	Vou[one,two]<-output[2,j]
	Vou[two,one]<-output[2,j]
}

return(Vou)
}

