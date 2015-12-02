#ifndef __WRAPVECTOR_H__
#define __WRAPVECTOR_H__

#include <RcppCommon.h>
#include "Utilities/Vector.h"

//We start with the forward declarations:
namespace Rcpp{
// conversion from C++ to R
template <> SEXP wrap(const Purple::Vector& vc);
}
#endif

