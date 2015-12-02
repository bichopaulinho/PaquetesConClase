
#include <Rcpp.h>
#include "ObjectiveFunction/DeJongFunction.h"
#include "OptimizationAlgorithm/Gradient.h"

#include <iostream>
#include <time.h>
#include <stdexcept>

using namespace Rcpp;
using namespace Purple;


//' PruebaGradientDescend
//'
//' @export
//' @useDynLib GradDesc
// [[Rcpp::export]]
NumericVector PruebaGradientDescend() {


    // De Jong function object

    DeJongFunction deJongFunction;

    int numberOfVariables = 3;

    deJongFunction.setNumberOfVariables(numberOfVariables);

    Vector<double> lowerBound(numberOfVariables, -5.12);
    Vector<double> upperBound(numberOfVariables, -5.12);

    deJongFunction.setLowerBound(lowerBound);
    deJongFunction.setUpperBound(upperBound);

    // Gradient descent object

    GradientDescent gradientDescent(&deJongFunction);

    gradientDescent.setOptimalStepSizeMethod(GradientDescent::BrentMethod);

    gradientDescent.setOptimalStepSizeTolerance(1.0e-6);

    gradientDescent.setEvaluationGoal(0.0);
    gradientDescent.setGradientNormGoal(0.0);
    gradientDescent.setMaximumNumberOfIterations(100);
    gradientDescent.setMaximumTime(1000.0);

    gradientDescent.setShowPeriod(10);

    Vector<double> initialArgument(numberOfVariables, 1.0);

    gradientDescent.setInitialArgument(initialArgument);

    gradientDescent.print();

    gradientDescent.save("GradientDescent.dat");

    Vector<double> minimalArgument = gradientDescent.getMinimalArgument();

    gradientDescent
        .saveOptimizationHistory("OptimizationHistory.dat");

    std::cout << std::endl;

    //return(0);
    return(Rcpp::wrap<Purple::Vector<double>>(minimalArgument));

}
