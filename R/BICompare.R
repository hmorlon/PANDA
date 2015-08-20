BICompare <- function(phylo,t,meth=c("ultrametric")){
	options(warn=-1)
	#get lLk, BIC
	kmeansBIC <- function(fit){
		m <- ncol(fit$centers)
		n <- length(fit$cluster)
		k <- nrow(fit$centers)
		D <- fit$tot.withinss
	return(data.frame(BIC = D + log(n)*m*k))
	}
		phyloM <- as.matrix(dist.nodes(phylo))
		rDP <- c()
		c()->r
		c()->p
		c()->q	
	
	##if tree is non-ultrametric	
		if(meth=="non-ultrametric"){	
			c() -> rDP
			c() -> r
			for(i in c(1:100)){
				rDP[[i]] <- as.matrix(
					dist.nodes(
						rtree(n=length(phylo$tip.label),
							br=runif(100,min=min(phylo$edge.length),
						max=max(phylo$edge.length)))))
					r[i] <- kmeansBIC(kmeans(rDP[[i]],t,algorithm="Hartigan-Wong"))
				}
			##get z-score for BICs
				c() -> zDP
				c() -> pDP
				m = mean(as.numeric(r))
				s = sd(as.numeric(r))
				for(j in 1:length(r)){
					pDP[[j]] <- 1-pnorm(as.numeric(r[j]),m,s)
				}
			col_r <- cbind(as.numeric(r),pDP)
		low_p <- subset(col_r,col_r[,2] == min(col_r[,2]))	
	r = low_p[,1]			
	}		

	##if tree is ultrametric	
		if(meth=="ultrametric"){
			c() -> rDP
			c() -> r
			for(i in c(1:100)){
				rDP[[i]] <- as.matrix(
					dist.nodes(
						rcoal(n=length(phylo$tip.label),
							br=runif(100,0,max(branching.times(phylo)/20)))))
			r[i] <- kmeansBIC(kmeans(rDP[[i]],t,algorithm="Hartigan-Wong"))
				}
			##get z-score for BICs
				c() -> zDP
				c() -> pDP
				m = mean(as.numeric(r))
				s = sd(as.numeric(r))
				for(j in 1:length(r)){
					pDP[[j]] <- 1-pnorm(as.numeric(r[j]),m,s)
				}
			col_r <- cbind(as.numeric(r),pDP)
		low_p <- subset(col_r,col_r[,2] == min(col_r[,2]))	
	r = low_p[,1]			
	}	
	
	##get BICs for tree and control
			kmeans(phyloM,t,algorithm="Hartigan-Wong")->q
				kmeansBIC(q)->p
rp<-cbind(p,r)
colnames(rp)<-c("tree BIC","random BIC")
res<-list("BIC_test"=rp,"clusters"=q$cluster,"BSS/TSS"=q$betweenss/q$totss)

class(res)	<- "BICompare"
return(res)

}