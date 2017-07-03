spectR<-function (phylo, method = c("standard")) 
{
    kurtosis.sub <- function(x, na.rm = FALSE, method = c("moment"), 
        ...) {
        method = match.arg(method)
        stopifnot(NCOL(x) == 1)
        if (!is.numeric(x) && !is.complex(x) && !is.logical(x)) {
            warning("argument is not numeric or logical: returning NA")
            return(as.numeric(NA))
        }
        if (na.rm) 
            x = x[!is.na(x)]
        n = length(x)
        if (is.integer(x)) 
            x = as.numeric(x)
        if (method == "moment") {
            kurtosis = sum((x - mean(x))^4/as.numeric(var(x))^2)/length(x)
        }
        if (method == "excess") {
            kurtosis = sum((x - mean(x))^4/var(x)^2)/length(x) - 
                3
        }
        if (method == "fisher") {
            kurtosis = ((n + 1) * (n - 1) * ((sum(x^4)/n)/(sum(x^2)/n)^2 - 
                (3 * (n - 1))/(n + 1)))/((n - 2) * (n - 3))
        }
        kurtosis
    }
    skewness <- function(x, na.rm = FALSE) {
        if (is.matrix(x)) 
            apply(x, 2, skewness, na.rm = na.rm)
        else if (is.vector(x)) {
            if (na.rm) 
                x <- x[!is.na(x)]
            n <- length(x)
            (sum((x - mean(x))^3)/n)/(sum((x - mean(x))^2)/n)^(3/2)
        }
        else if (is.data.frame(x)) 
            sapply(x, skewness, na.rm = na.rm)
        else skewness(as.vector(x), na.rm = na.rm)
    }
        integr <- function(x, f) {
            if (!is.numeric(x)) {
                stop("The variable of integration \"x\" is not numeric.")
            }
            if (!is.numeric(f)) {
                stop("The integrand \"f\" is not numeric.")
            }
            if (length(x) != length(f)) {
                stop("The lengths of the variable of integration and the integrand do not match.")
            }
            n = length(x)
            integral = 0.5 * sum((x[2:n] - x[1:(n - 1)]) * (f[2:n] + 
                f[1:(n - 1)]))
            return(integral)
        }
        
    
    if (method == "standard") {
        e = eigen(graph.laplacian(graph.adjacency(data.matrix(dist.nodes(phylo)), 
            weighted = T), normalized = F), symmetric = T, only.values = F)
        m = subset(e$values, e$values >= 1)
        d = dens(log(m))
        dsc = d$y/integr(d$x, d$y)
        gaps <- abs(diff(m))
        gapMat <- as.matrix(gaps)
        modalities <- c(1:length(gapMat))
        gapMatCol <- cbind(modalities, gapMat)
        eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
            2]))
        principal_eigenvalue <- max(m)
        kurtosis <- kurtosis.sub(dsc)
        skewness <- skewness(dsc)
        peak_height <- max(dsc)
        res <- list(eigenvalues = e$values, principal_eigenvalue = principal_eigenvalue, 
            asymmetry = skewness, peakedness1 = kurtosis, peakedness2 = peak_height, 
            eigengap = eigenGap[, 1])
    }
    if (method == "normal") {
        e = eigen(graph.laplacian(graph.adjacency(data.matrix(dist.nodes(phylo)), 
            weighted = T), normalized = T), symmetric = T, only.values = F)
        m = subset(e$values, e$values >= 0)
        d = dens(m)
        dsc = d$y/integr(d$x, d$y)
        gaps <- abs(diff(m))
        gapMat <- as.matrix(gaps)
        modalities <- c(1:length(gapMat))
        gapMatCol <- cbind(modalities, gapMat)
        eigenGap <- subset(gapMatCol, gapMatCol[, 2] == max(gapMatCol[, 
            2]))
        principal_eigenvalue <- max(m)
        kurtosis <- kurtosis.sub(dsc)
        skewness <- skewness(dsc)
        peak_height <- max(dsc)
        res <- list(eigenvalues = e$values, principal_eigenvalue = principal_eigenvalue, 
            asymmetry = skewness, peakedness1 = kurtosis, peakedness2 = peak_height, 
            eigengap = eigenGap[, 1])
    }
    class(res) <- "spectR"
    return(res)
}
