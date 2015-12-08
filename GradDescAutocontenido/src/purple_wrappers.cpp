#include "purple_wrappers.h"
#include <Rcpp.h>
using namespace Rcpp;

namespace Rcpp{

    template <> SEXP wrap(const Purple::Vector<double>& vc){
        //Especialización de la plantilla de wrap para transformar Purple::Vector en Rcpp::NumericVector
        NumericVector res(vc.begin(),vc.end());
        return res;
    }
    template <> Purple::Vector<double> as(SEXP vector_r){
        Rcpp::NumericVector vect(vector_r);
        Purple::Vector<double> p_vect(vect.size());
        std::copy(vect.begin(),vect.end(),p_vect.begin());
        return p_vect;
    }

    template <> Purple::Matrix<double> as(SEXP matrix_r){
        Rcpp::NumericMatrix A(matrix_r);
        int ncol = A.ncol();
        int nrow = A.nrow();

        Purple::Matrix<double> M_A(nrow,ncol);
        for (int i = 0;i<nrow;i++){
            Rcpp::NumericMatrix::Row Arow = A(i, _);
            std::copy(Arow.begin(),Arow.end(),M_A[i]);
        }
        return M_A;
    }
}
