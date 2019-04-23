#  This class represents a matrix which is the product of a diagonal matrix
#  by a Toeplitz matrix: A = D0 %*% Toeplitz.
#  This special structure allows fast matrix-vector and exponential computations.
#
#  This class defines generic methods applyV and expV, and concrete implementations
#  are defined in derived classes MPhiFFT and MPhiNoFFT
setClass(
    Class = "MPhiAbstract"
)

setGeneric(
    name="applyV",
    def=function(object, x){standardGeneric("applyV")}
)

setGeneric(
    name="expV",
    def=function(object, diag1, diag2, x){standardGeneric("expV")}
)

