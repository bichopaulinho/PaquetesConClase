#include "Utilities/Vector.h"
#include "wrap_purple_vector.h"
#include <Rcpp.h>
using namespace Rcpp;

namespace Rcpp{

    template <> SEXP wrap(const Purple::Vector<double>& vc){
        //Especialización de la plantilla de wrap para transformar Purple::Vector en Rcpp::NumericVector
        NumericVector res(vc.getSize());

        for (auto i = vc.begin(); i != vc.end(); i++){
            res[i] = vc[i];
        }
        return res;
    }

}
