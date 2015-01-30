list1<-paste("d",rep(LETTERS[1:length(phylo$tip.label)],each=length(phylo$tip.label)),".",seq(1,length(phylo$tip.label)),sep="") #this names the variables for the n+n^2 series of ODEs needed, WILL NEED TO MODIFY FOR MORE THAN 26 tips
C<-matrix(list1,nrow=length(phylo$tip.label)) #coerces this list into a matrix
cov.list<-C[lower.tri(C,diag=FALSE)]#this extracts covariance terms
var.list<-diag(C) #this names variance terms separately

##define ODEs for deSOLVE
 ou<-function(t, state, parameters) {
 with(as.list(c(state, parameters)),{
 # rate of change

	##this loop defines the set of ODEs for the variance terms (stored in "var.list"), one for each tip of the phylo
	##note: if you want to see visualize the output of this, replace the line beginning with 'eval(' with:
	#assign(var.list[i],paste("(-(((2*(1-length(phylo$tip.label)))/length(phylo$tip.label)*s)+2*a)*",float,")+(((2*s)/length(phylo$tip.label))*sum(",sumterm,"))+b",sep=""))
	##then call each term in var.list (e.g., dA.1) to see what the ODE for that term is
	for(i in 1:length(var.list)){
		float<-substr(var.list[i], 2, 12) #this removes the d
		##need to coerce sumterm into lower triangle
		sumint<-vector() 
		for(j in 1:length(setdiff(1:length(phylo$tip.label),match(substr(var.list[i], 2,2),LETTERS)))){
			if(setdiff(1:length(phylo$tip.label),match(substr(var.list[i], 2,2),LETTERS))[j]>match(substr(var.list[i], 2,2),LETTERS)){
			sumint[j]<-C[setdiff(1:length(phylo$tip.label),match(substr(var.list[i], 2,2),LETTERS))[j],match(substr(var.list[i], 2,2),LETTERS)]} else {
			sumint[j]<-C[match(substr(var.list[i], 2,2),LETTERS),setdiff(1:length(phylo$tip.label),match(substr(var.list[i], 2,2),LETTERS))[j]]
			}
			}
		sumterm<-paste(substr(sumint,2,12),collapse=",")		
	 	eval(parse(text=paste(var.list[i],"<-(-(((2*(1-length(phylo$tip.label)))/length(phylo$tip.label)*s)+2*a)*",float,")+(((2*s)/length(phylo$tip.label))*sum(",sumterm,"))+b",sep="")))
	 	}	
	
	##this loop defines the set of ODEs for the covariance terms (stored in "cov.list"), only for the lower triangle elements (including diagonal) of the VCV matrix
	##note: if you want to see visualize the output of this, replace the line beginning with 'eval(' with:
	#assign(cov.list[i],paste("(-(((2*(1-length(phylo$tip.label)))/length(phylo$tip.label)*s)+2*a)*",float,")+(((s)/length(phylo$tip.label))*sum(",sumterm,"))",sep=""))
	##then call each term (e.g., dA.2) to see what the ODE for that term is
	for(i in 1:length(cov.list)){
		float<-substr(cov.list[i], 2, 12) #this removes the d
		sumterms.raw<-setdiff(C[,c(as.numeric(match(substr(cov.list[i], 2,2),LETTERS)),as.numeric(substr(cov.list[i], 4,12)))],c(C[,c(as.numeric(match(substr(cov.list[i], 2,2),LETTERS)),as.numeric(substr(cov.list[i], 4,12)))][as.numeric(substr(cov.list[i], 4,12)),1],C[,c(as.numeric(match(substr(cov.list[i], 2,2),LETTERS)),as.numeric(substr(cov.list[i], 4,12)))][as.numeric(match(substr(cov.list[i], 2,2),LETTERS)),2]))
		sumint<-vector()
		for(j in 1:length(sumterms.raw)){
			if((sumterms.raw[j] %in% cov.list)||(sumterms.raw[j] %in% var.list)){
			sumint[j]<-sumterms.raw[j]} else { 
			sumint[j]<-C[ceiling(match(sumterms.raw[j],C)/length(phylo$tip.label)),match(sumterms.raw[j],C)%%length(phylo$tip.label)]			
			}
			}
			if((if(match(cov.list[i],C)%%length(phylo$tip.label)==0){length(phylo$tip.label)} else {match(cov.list[i],C)%%length(phylo$tip.label)})==ceiling(match(cov.list[i],C)/length(phylo$tip.label))){
			sumterm<-paste(substr(sumint,2,12),substr(sumint,2,12),sep=",",collapse=",")} else {sumterm<-paste(substr(sumint,2,12),collapse=",")}				
			eval(parse(text=paste(cov.list[i],"<-(-(((2*(1-length(phylo$tip.label)))/length(phylo$tip.label)*s)+2*a)*",float,")+(((s)/length(phylo$tip.label))*sum(",sumterm,"))",sep="")))
		}

				
 # return the rate of change
 sttrt=paste(c(var.list,cov.list),sep="",collapse=",")
 sttat=paste("list(c(",sttrt,"))",sep="")
 eval(parse(text=sttat))
 }) # end with(as.list ...
 }
 
parameters<-c(a=alpha,b=sigma,s=sterm) 

##the next three lines give starting values (all 0 here, as indicated by "=0", can modify in the future, perhaps using
##x<-match("C.3",row.names(as.data.frame(state)))
##state[x]<-new.value

start=c(paste(substr(c(var.list,cov.list),2,14),"=0",sep="",collapse=","))
state.new<-paste("c(",start,")",sep="")
state<-eval(parse(text=state.new))


U<-unique(c(vcv.phylo(phylo)))
output<-ode(y=state,times=sort(U),func=ou,parms=parameters)

