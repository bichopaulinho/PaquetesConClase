#include "wrap_purple_vector.h"
#include <Rcpp.h>
#include "ObjectiveFunction/SumaCuadrados.h"
#include "OptimizationAlgorithm/GradientDescent.h"

#include <iostream>
#include <time.h>
#include <stdexcept>

using namespace Rcpp;
//using namespace Purple;


//' PruebaGradientDescent
//'
//' @export
//' @useDynLib GradDesc
// [[Rcpp::export]]
Rcpp::NumericVector GradientDescentSumaCuadrados(Rcpp::NumericMatrix & A,const Rcpp::NumericVector & b) {

    int ncol = A.ncol();
    int nrow = A.nrow();
    if (b.size() != nrow){
        Rcpp::stop("El vector b tiene que tener la misma longitud que las filas de A");
    }

    // Copiamos la matriz A
    Purple::Matrix<double> M_A(nrow,ncol);
    for (int i = 0;i<nrow;i++){
        Rcpp::NumericMatrix::Row Arow = A(i, _);
        std::copy(Arow.begin(),Arow.end(),M_A[i]);
    }

    // Copiamos el vector b
    Purple::Vector<double> M_b(b.size());
    std::copy(b.begin(),b.end(),M_b.begin());

    Purple::SumaCuadrados funcion_objetivo(M_b,M_A);
    // Gradient descent object

    Purple::GradientDescent gradientDescent(&funcion_objetivo);

    gradientDescent.setOptimalStepSizeMethod(Purple::GradientDescent::BrentMethod);

    gradientDescent.setOptimalStepSizeTolerance(1.0e-6);

    gradientDescent.setEvaluationGoal(0.0);
    gradientDescent.setGradientNormGoal(1.0e-6);
    gradientDescent.setMaximumNumberOfIterations(100);
    gradientDescent.setMaximumTime(1000.0);

    gradientDescent.setShowPeriod(10);

    Purple::Vector<double> initialArgument(ncol, 1.0);

    gradientDescent.setInitialArgument(initialArgument);

    gradientDescent.print();

    //gradientDescent.save("GradientDescent.dat");

    Purple::Vector<double> minimalArgument = gradientDescent.getMinimalArgument();

    //gradientDescent
    //    .saveOptimizationHistory("OptimizationHistory.dat");


    return(Rcpp::wrap<Purple::Vector<double>>(minimalArgument));
}
