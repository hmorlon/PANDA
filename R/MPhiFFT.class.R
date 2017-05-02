library(pracma)

#  This class represents a matrix A = (1/rowSums(Toep)) * Toep
#  where Toep is a Toeplitz matrix.
#  This special structure allows fast matrix-vector and exponential computations.
#
#  A Toeplitz matrix T of size NxN can be embedded into a circulant matrix of
#  size 2Nx2N  [ T S ]
#              [ S T ]
#  This circulant matrix is fully determined by its first column c.
#  Then T*x = (I  0) A [ x ]
#                      [ 0 ]
#  and A*X0 = IFFT(FFT(c) * FFT(X0))
#  In order to improve FFT performance, we enlarge A such that its dimension
#  is a power of two.

MPhiFFT <- setClass(
    "MPhiFFT",
    slots = c(
        method= "character",
        d0 = "vector",
        firstColumn = "vector",
        firstRow = "vector",
        fftColumn = "vector",
        n = "numeric",
        zero = "numeric"
    ),
    contains="MPhiAbstract"
)


setMethod(
    f="initialize",
    signature=c("MPhiFFT"),
    definition=function(.Object, firstColumn, firstRow, ...){
      callNextMethod()
      args <- list(...)
      if(is.null(args$method)) .Object@method="FFT"
      else .Object@method=args$method
      .Object@firstColumn <- firstColumn
      .Object@firstRow <- firstRow
      if(length(.Object@firstColumn) != length(.Object@firstRow)) {
          stop("Error: first column and first row must have the same size")
      }
      if(.Object@firstColumn[1] != .Object@firstRow[1]) {
          stop("Error: first row and first column must have the same first element")
      }
      .Object@n <- length(firstColumn)
      .Object@d0 <- sapply(1:length(firstColumn), function(irow) {return(1/sum(c(firstColumn[irow:2], firstRow[1:(.Object@n+1-irow)])))})
      # Compute the next higher power of 2 for FFT, it is better to perform FFT on
      # vectors of length power of 2.
      N2 = 2 * .Object@n
      size=2^ceiling(log(N2)/log(2))
      circulantColumn=c(firstColumn, rep(0, size + 1 - 2*.Object@n), firstRow[length(firstRow):2])
      .Object@zero <- rep(0, size - .Object@n)
      .Object@fftColumn <- Re(fft(circulantColumn))
      return(.Object)
    }
)

setMethod(
    f="applyV",
    signature=c("MPhiFFT","vector"),
    definition=function(object,x){
      return(object@d0 * Re(pracma::ifft(object@fftColumn * fft(c(x,object@zero)))[1:object@n]))
    }
)

# Compute exp(diag(diag1)+M*diag2)%*%x = exp(diag(diag1)+diag(diag2)%*%M)%*%x
#                                      = exp(diag1*x + (diag2*d0)*(toeplitz%*%x))
setMethod(
    f="expV",
    signature=c("MPhiFFT","vector","vector","vector"),
    definition=function(object,diag1,diag2,x){
      mask_diag1 = diag1 < -20
      x[mask_diag1] = 0
      expGv=x
      Gnv=x
      i=1
      epsnormv = 1e-10 * sum(abs(x))
      MATVECT <- selectMethod(applyV, c("MPhiFFT","vector"))
      while (sum(abs(Gnv)) > epsnormv){
        Gnv = (diag1*Gnv + diag2*MATVECT(object,Gnv)) / i
        Gnv[mask_diag1] = 0
        expGv = expGv + Gnv
        i=i+1
      }
      expGv[mask_diag1] = 2e-9
      return(expGv)
    }
)

setMethod(
    f="[",
    signature=c("MPhiFFT","numeric","missing","missing"),
    definition=function(x,i,j,..., drop=FALSE){
      if (nargs() == 3) return(x@d0[i]*c(x@firstColumn[i:2], x@firstRow[1:(x@n+1-i)]))
      else stop("MPhiFFT[i] not implemented")
    }
)

