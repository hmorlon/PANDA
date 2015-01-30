#include <Rcpp.h>
#include <iostream>
#include <vector>

using namespace Rcpp;

// [[Rcpp::export]]
void fillIndices(IntegerVector & csrJ, IntegerVector & csrP, NumericVector & csrX, int len) {

  int pos = 0;
  int dim = len*(len+1)/2;
  std::vector<int> indI(dim);
  std::vector<int> indJ(dim);
  for (int k = 0; k < len; ++k)
  {
    indI[k] = k;
    indJ[k] = k;
  }
  int counter = len;
  for (int l = 0; l < len-1; ++l)
      for (int k = l+1; k < len; ++k, ++counter)
      {
        indI[counter] = k;
        indJ[counter] = l;
      }
  for (int ij = 0; ij < dim; ++ij)
  {
    const int i0(indI[ij]);
    const int j0(indJ[ij]);
    for (int kl = 0; kl < dim; ++kl)
    {
      const int k0(indI[kl]);
      const int l0(indJ[kl]);
      if ((ij != kl) && (k0 == i0 || k0 == j0 || l0 == i0 || l0 == j0)) {
        csrJ[pos++] = kl + 1;
      }
    }
  }
  for (int i = 0; i <= len; ++i)
    csrP[i] = i * (len - 1);
  for (int i = len + 1; i <= dim; ++i)
    csrP[i] = (2 * i - len) * (len - 1);
  std::fill(csrX.begin(), csrX.begin() + len*(len-1), 2.0);
  std::fill(csrX.begin() + len*(len-1), csrX.end(), 1.0);
}

