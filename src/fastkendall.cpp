#include <algorithm>
#include <vector>
#include <stdint.h>

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>

#include "goodman_kruskal_gamma.hpp"


template<typename compare_t = std::less<double>, typename T = std::vector<double> >
struct compare_other{
  // Functor for use with std::sort.
  // Don't compare directly but use i, j as indices for another vector/array
  compare_other(const T& other_): other(&other_), compare(compare_t()) {}

  bool operator()(size_t i, size_t j) const{
    return compare((*other)[i], (*other)[j]);
  }

private:
  const T* other;
  compare_t compare;
};  


extern "C"{
  SEXP concordance_count(SEXP X, SEXP Y){
    // prepare for and do concordance_count() from goodman_kruskal_gamma.hpp,
    // returns numeric vector
    // c(concordant, discordant, extraX, extraY, spare)

    size_t N = length(X);	// test if length(X) == length(Y) etc. in R
    double* Xptr = REAL(X);
    double* Yptr = REAL(Y);

    // for(size_t i = 0; i != N; ++i)
    //   std::cout << Xptr[i] << " ";
    // std::cout << "\n";
    //   for(size_t i = 0; i != N; ++i)
    // 	std::cout << Yptr[i] << " ";
    //   std::cout << "\n";
    
    // preparation
    // sort X and Y
    std::vector<size_t> order_X(N); // use uint32_t if memory becomes a problem
    // get order with compare_other
    for(size_t i = 0; i != order_X.size(); ++i)
      order_X[i] = i;
    sort(order_X.begin(), order_X.end(), compare_other<std::less<double>, double*>(Xptr));

    // reorder X and Y
    std::vector<double> X_ord(N);
    std::vector<double> Y_ord(N);
    for(size_t i = 0; i != N; ++i)
      X_ord[i] = Xptr[order_X[i]];
    for(size_t i = 0; i != N; ++i)
      Y_ord[i] = Yptr[order_X[i]];


    // sort Y for equal values of X
    secondary_sort(X_ord.begin(), X_ord.end(), Y_ord.begin(), Y_ord.end());
    
    
    // reserve result
    SEXP result;
    PROTECT( result = allocVector(REALSXP, 5) );
    uint64_t tmp[5];
    
    // calculate result
    concordance_count(X_ord.begin(), X_ord.end(), Y_ord.begin(), Y_ord.end(),
		      tmp[0], tmp[1], tmp[2], tmp[3], tmp[4]);

    std::copy(tmp, tmp+5, REAL(result));
    
    UNPROTECT(1);
    return result;
  }
}
