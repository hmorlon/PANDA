#plot spectral density
spectR <- function(phylo,method=c("standard")){

##kurtosis
kurtosis.sub <-
    function (x, na.rm = FALSE, method = c("moment"), ...)
{
    
    method = match.arg(method)

    stopifnot(NCOL(x) == 1)

    # Warnings:
    if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
        warning("argument is not numeric or logical: returning NA")
        return(as.numeric(NA))}

    # Remove NAs:
    if (na.rm) x = x[!is.na(x)]

    # Kurtosis:
    n = length(x)
    if (is.integer(x)) x = as.numeric(x)
    if (method == "moment") {
        kurtosis = sum((x-mean(x))^4/as.numeric(var(x))^2)/length(x)
    }
     if (method == "excess") {
        kurtosis = sum((x-mean(x))^4/var(x)^2)/length(x) - 3
    }

    if (method == "fisher") {
        kurtosis = ((n+1)*(n-1)*((sum(x^4)/n)/(sum(x^2)/n)^2 -
            (3*(n-1))/(n+1)))/((n-2)*(n-3))
    }

    # Return Value:
    kurtosis
}


#skewness
skewness <- function (x, na.rm = FALSE) 
{
    if (is.matrix(x)) 
        apply(x, 2, skewness, na.rm = na.rm)
    else if (is.vector(x)) {
	if (na.rm) x <- x[!is.na(x)] 
	n <- length(x)
     (sum((x-mean(x))^3)/n)/(sum((x-mean(x))^2)/n)^(3/2)
	}
    else if (is.data.frame(x)) 
        sapply(x, skewness, na.rm = na.rm)
    else skewness(as.vector(x), na.rm = na.rm)
}
		
	if(method=="standard"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=F)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=1)

	#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
		
		#get principal eigenvalue
		max(m) -> principal_eigenvalue
					
		#get kurtosis
		kurtosis.sub(m) -> kurtosis
		
		#get skewness
		skewness(m) -> skewness
		
		#output
		res<-list(eigenvalues=e$values, principal_eigenvalue=principal_eigenvalue,asymmetry=skewness,peakedness=kurtosis, eigengap=eigenGap[,1])
	}
	
	if(method=="normal1"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=T)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=0)
	

		#get eigengap
		abs(diff(m))->gaps
			as.matrix(gaps)->gapMat
				c(1:length(gapMat))->modalities
			cbind(modalities,gapMat)->gapMatCol
		subset(gapMatCol,gapMatCol[,2]==max(gapMatCol[,2]))->eigenGap
				
		#get principal eigenvalue
		max(m) -> principal_eigenvalue
					
		#get kurtosis
		kurtosis.sub(m) -> kurtosis
		
		#get skewness
		skewness(m) -> skewness
		
		#output
		res<-list(eigenvalues=e$values, principal_eigenvalue=principal_eigenvalue,asymmetry=skewness,peakedness=kurtosis,eigengap=eigenGap[,1])			
			}

	if(method=="normal2"){
		e=eigen(
			graph.laplacian(
				graph.adjacency(
					data.matrix(
						dist.nodes(phylo))
					,weighted=T)
				,normalized=F)
			,symmetric=T,only.values=F)
		m=subset(e$values,e$values>=0)

		#get principal eigenvalue
		max(m) -> principal_eigenvalue
		
		#get skewness
		skewness(m) -> skewness
					
		#get kurtosis
		kurtosis.sub(m) -> kurtosis
		
				#output
		res<-list(eigenvalues=e$values,principal_eigenvalue=principal_eigenvalue,asymmetry=skewness,peakedness=kurtosis,eigengap=eigenGap[,1])

	}
	class(res)	<- "spectR"
	return(res)					
}


