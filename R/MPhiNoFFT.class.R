library(expoRkit)
library(expm)
library(Matrix)

#  This class represents a matrix B = (1/rowSums(Toep)) * Toep
#  where Toep is a Toeplitz matrix.
MPhiNoFFT <- setClass(
    Class = "MPhiNoFFT",
    slots = c(
        B= "matrix",                  # Dense representation
        sparse= "logical",            # Whether dense representation is stored in sparse format or not
        method= "character",          # Method to use in expV
        applyV = "function",          # Compute A %*% x
        expV = "function"             # Compute exp(diag1 + A*diag2) %*% x
    ),
    contains="MPhiAbstract"
)


setMethod(
    f="initialize",
    signature=c("MPhiNoFFT"),
    definition=function(.Object, B, ...){
      callNextMethod()
      .Object@B <- B
      .Object@sparse <- FALSE
      args <- list(...)
      if(is.null(args$method)) .Object@method="Rexpv"
      else .Object@method=args$method

      return(.Object)
    }
)

setMethod(
    f="applyV",
    signature=c("MPhiNoFFT","vector"),
    definition=function(object,x){
      return(object@B %*% x)
    }
)

setMethod(
    f="expV",
    signature=c("MPhiNoFFT","vector","vector","vector"),
    definition=function(object,diag1,diag2,x){
      B = diag(diag1) + object@B * diag2
      if (object@sparse) {B = Matrix(B, sparse=T)}
      switch( EXPR=object@method,
              "Sidje98"={return(expm::expAtv(B,x)$eAtv)},
              "expoRkit"={return(expoRkit::expv(x=B,v=x))},
              "Rexpv"={return(expoRkit::Rexpv(B@x,B@i+1,B@p+1,length(x),x))},
              stop("Invalid method name") )
    }
)

setMethod(
    f="[",
    signature=c("MPhiNoFFT","numeric","missing","missing"),
    definition=function(x,i,j,..., drop=FALSE){
      if (nargs() == 3) return(x@B[i,])
      else stop("MPhiNoFFT[i] not implemented")
    }
)

