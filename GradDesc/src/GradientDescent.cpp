#include "purple_wrappers.h"
#include <Rcpp.h>
#include "ObjectiveFunction/SumaCuadrados.h"
#include "OptimizationAlgorithm/GradientDescent.h"

#include <iostream>
#include <time.h>
#include <stdexcept>

using namespace Rcpp;


//' GradientDescentSumaCuadrados
//'
//' @param A Matriz con regresores
//' @param b Vector objetivo
//'
//' @detail
//'
//' @export
//' @useDynLib GradDesc
// [[Rcpp::export]]
Rcpp::NumericVector GradientDescentSumaCuadrados(SEXP A,SEXP b) {

    // Copiamos el vector b
    Purple::Vector<double> M_b = Rcpp::as<Purple::Vector<double>>(b);
    Purple::Matrix<double> M_A = Rcpp::as<Purple::Matrix<double>>(A);
    int ncol = M_A.getNumberOfColumns();
    int nrow = M_A.getNumberOfRows();
    if (M_b.getSize() != nrow){
        Rcpp::stop("El vector b tiene que tener la misma longitud que las filas de A");
    }


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
