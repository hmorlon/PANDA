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
##nodeDist, the distance between each node and the root
##nodeDiff, the distance between each node and the previous node (times for integration)

nodeDist<-vector()
root <- length(phylo$tip.label) + 1
heights<-nodeHeights(phylo)
for (i in root:(length(phylo$tip.label)+phylo$Nnode)){
	if(i==root){
	nodeDist[i-length(phylo$tip.label)]<-0
	} else {
	int<-heights[match(i,phylo$edge[,1]),1]
	nodeDist[i-length(phylo$tip.label)]<-int
	}
}
nodeDist<-c(nodeDist,max(heights))
nodeDiff<-diff(nodeDist)
##label the branches for each segment of tree to be integrated and identify the node at which the branch terminates

mat<-matrix(ncol=3)
for(i in 1:phylo$Nnode){
	branches<-phylo$edge[phylo$edge[,1]==i+length(phylo$tip.label),]
	for(j in 1:2){
		if(branches[j,2]>length(phylo$tip.label)){
			int<-matrix(ncol=3)
			int[1]<-i+length(phylo$tip.label)
			int[2]<-paste(".",THREELETTERS[1],sep="")
			int[3]<-phylo$edge[which(phylo$edge[,1]==i+length(phylo$tip.label)),][j,2]
			k<-1
			while(!is.na(match(int[2],mat))){
				int[2]<-paste(".",THREELETTERS[k],sep="")
				k<-k+1
				}
			} else {
			int<-matrix(ncol=3)
			int[1]<-i+length(phylo$tip.label)
			int[2]<-phylo$tip.label[branches[j,2]]
			int[3]<-0 ##NOTE :I am considering tips to be "0" here and use this below
			}
		if(is.na(mat[1])){
		mat<-int} else {
		mat<-rbind(mat,int)}
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

for(i in 1:phylo$Nnode){ ##for each interval between branches...

##var.list is the list of terms for which there is a variance term at each step in time
##cov.list is the list of covariance terms
var.list<-unlist(nat[[i]])
cov.list<-apply(t(combn(unlist(nat[[i]]),2)),1,paste,collapse="_")
Cmat<-matrix(nrow=length(var.list),ncol=length(var.list)) ##"Cmat" is a lower triangle matrix with each of the variance and covariance terms, used to define equations below
diag(Cmat)<-var.list
Cmat[lower.tri(Cmat)]<-cov.list

##define ODEs for deSOLVE
 ou<-function(t, state, parameters) {
 with(as.list(c(state, parameters)),{
 # rate of change

	##this loop defines the set of ODEs for the variance terms (stored in "var.list"), one for each "extant branch" of the phylo
	##note: if you want to visualize the output of this, replace the line beginning with 'eval(' with:
	##then call each term in var.list (i.e., in var.list) to see what the ODE for that term is
	for(m in 1:length(var.list)){
		sumint<-intersect(Cmat[,m],Cmat[lower.tri(Cmat)])
		if(m>1){ #unless it is the first term in varlist, it is necessary to coerce the script to extract the lower triangle version of a covariance term (e.g., AB instead of BA)
			sumint<-c(sumint, Cmat[m,1:(m-1)])
			}
		sumterm<-paste(sumint,collapse=",")		
	 	eval(parse(text=paste("d",var.list[m],"<--(((((2*(length(var.list)-1))/length(var.list)*s))+(2*a))*",var.list[m],")+(((2*s)/length(var.list))*sum(",sumterm,"))+b",sep="")))
	 	#assign(var.list[m],paste("-(((((2*(length(var.list)-1))/length(var.list)*s))+(2*a))*",var.list[m],")+(((2*s)/length(var.list))*sum(",sumterm,"))+b",sep=""))
	 	}	
	 

	##this loop defines the set of ODEs for the covariance terms (stored in "cov.list"), only for the lower triangle elements (including diagonal) of the VCV matrix
	##elements included for each covariance include ALL variance and covariance terms (other than the covariance being computed) involving the two component terms
	##note: if you want to visualize the output of this, replace the line beginning with 'eval(' with:
	##then call each term (i.e., in cov.list) to see what the ODE for that term is
	for(m in 1:length(cov.list)){
		sumint<-vector()
		for(j in 1:2){
		k<-which(Cmat==(unlist(strsplit(cov.list[m],"_"))[j]),arr.in=TRUE)[1]
		sumint<-c(sumint,intersect(Cmat[,k],Cmat[lower.tri(Cmat,diag=TRUE)]))
		if(k>1){
			sumint<-c(sumint, Cmat[k,1:(k-1)])
			}
		}
		sumint<-unique(sumint)[which(unique(sumint)!=cov.list[m])]
		sumterm<-paste(sumint,collapse=",")	##NEED TO FIND OUT if the order of terms is important in these sum terms, if not no worries IF SO then this needs some work
		eval(parse(text=paste("d",cov.list[m],"<-(-(((2*(length(var.list)-1))/length(var.list)*s)+2*a)*",cov.list[m],")+(((s)/length(var.list))*sum(",sumterm,"))",sep="")))
		#assign(cov.list[m],paste("(-(((2*(length(var.list)-1))/length(var.list)*s)+2*a)*",cov.list[m],")+(((s)/length(var.list))*sum(",sumterm,"))",sep=""))
		}
				
 # return the rate of change
 sttrt=paste("d",c(var.list,cov.list),sep="",collapse=",")
 sttat=paste("list(c(",sttrt,"))",sep="")
 eval(parse(text=sttat))
 }) # end with(as.list ...
 }
 
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
if(termlist[l] %in% colnames(output[[i-1]])){ #if the term to be numerically integrated was present in the previous generation(branch) of the tree
	tmin1<-output[[i-1]][2,match(termlist[l],colnames(output[[i-1]]))] #simply use the last value from that numerical integration as the starting value
	x<-match(termlist[l],row.names(as.data.frame(state)))  #and append that to the state vector
	state[x]<-tmin1
	}  
else { #if term is not present in previous generation/branch of the tree, lookup which branch it descends from
	if(termlist[l] %in% var.list){ ##this loop looks up values for variance terms using the "mat" matrix that is produced earlier
		prev.branch<-mat[mat[,3]==mat[mat[,2]==termlist[l],1],2]
		tmin1<-output[[i-1]][2,match(prev.branch,colnames(output[[i-1]]))]
		x<-match(termlist[l],row.names(as.data.frame(state))) 
		state[x]<-tmin1
	}
	else{ ##for covariance terms that weren't present in previous generation/branch, it is necessary to decompose covariance term, look up terms appropriately, and re-assemble them in the correct order so they match the lower triangle format used in definition of terms
		firstterm<-strsplit(termlist[l],"_")[[1]][1]
		scondterm<-strsplit(termlist[l],"_")[[1]][2]
		if(firstterm %in% colnames(output[[i-1]]) || scondterm %in% colnames(output[[i-1]])){
			if(firstterm %in% colnames(output[[i-1]])){
				prev.branch<-mat[mat[,3]==mat[mat[,2]==scondterm,1],2]
				prev.branch2<-firstterm} 
			else{  ##NOTE prev.branch2 just identifies which term was present previously
				prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
				prev.branch2<-scondterm
				}
		prev.term<-paste(prev.branch2,prev.branch,sep="_")
		if(is.na(match(prev.term,colnames(output[[i-1]])))){
			prev.term<-paste(prev.branch,prev.branch2,sep="_")
			}
		tmin1<-output[[i-1]][2,match(prev.term,colnames(output[[i-1]]))]
		x<-match(termlist[l],row.names(as.data.frame(state))) 
		state[x]<-tmin1
		} 				
	else {
		##neither first nor second term was present in the previous generation/branch, so the starting value is a variance term from t-1
		prev.branch<-mat[mat[,3]==mat[mat[,2]==firstterm,1],2]
		tmin1<-output[[i-1]][2,match(prev.branch,colnames(output[[i-1]]))]
		x<-match(termlist[l],row.names(as.data.frame(state))) 
		state[x]<-tmin1
	}
	}
	}
	}
	}


##NOW, run numerical integration for given time step and append to a list
output[[i]]<-ode(y=state,times=c(0,nodeDiff[i]),func=ou,parms=parameters)
}


#########	RETURN VCV MATRIX  	#########

Vou<-mrca(phylo)
diag(Vou)<-output[[length(nodeDiff)]][2,][2:(length(phylo$tip.label)+1)]
for(j in (length(phylo$tip.label)+2):length(output[[length(nodeDiff)]][2,])){

	string<-strsplit(names(output[[length(nodeDiff)]][2,j]),"_")
	string<-unlist(string)
	one<-match(string[1],colnames(Vou))
	two<-match(string[2],colnames(Vou))
	if(lower.tri(Vou)[one,two]==TRUE){
		Vou[one,two]<-output[[length(nodeDiff)]][2,j]} else{
		Vou[two,one]<-output[[length(nodeDiff)]][2,j]
	}
}

Vou[upper.tri(Vou)]<-t(Vou)[upper.tri(t(Vou))]
return(Vou)
}

