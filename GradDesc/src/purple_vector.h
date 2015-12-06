#ifndef __RCPP_PURPLE_VECTOR_H__
#define __RCPP_PURPLE_VECTOR_H__

#include <RcppCommon.h>
#include "Utilities/Vector.h"

//We start with the forward declarations:
// http://dirk.eddelbuettel.com/code/rcpp/Rcpp-extending.pdf
namespace Rcpp{
// conversion from C++ to R
    template <> SEXP wrap(const Purple::Vector<double>& vc);
    template <> Purple::Vector<double> as(SEXP);
}
#endif

