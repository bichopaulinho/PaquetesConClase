#include "Rcpp_purple_wrappers.h"
#include <Rcpp.h>
#include "ObjectiveFunction/SumaCuadrados.h"
#include "OptimizationAlgorithm/GradientDescent.h"


class GradientDescentSumaCuadradosClase {
public:
    GradientDescentSumaCuadradosClase(SEXP A, SEXP b):
                    funcion_objetivo(Rcpp::as<Purple::Vector<double>>(b),
                                     Rcpp::as<Purple::Matrix<double>>(A)),gradient_descent(){

        gradient_descent.setObjectiveFunction(&funcion_objetivo);
        gradient_descent.setOptimalStepSizeMethod(Purple::GradientDescent::BrentMethod);

        gradient_descent.setOptimalStepSizeTolerance(1.0e-6);

        gradient_descent.setEvaluationGoal(0.0);
        gradient_descent.setGradientNormGoal(1.0e-6);
        gradient_descent.setMaximumNumberOfIterations(1000);
        gradient_descent.setMaximumTime(1000.0);

        gradient_descent.setShowPeriod(10);

        Purple::Vector<double> initialArgument(funcion_objetivo.getNumberOfVariables(), 1.0);

        gradient_descent.setInitialArgument(initialArgument);

        gradient_descent.print();

    }
    void setEvaluationGoal(double ev_goal){
        gradient_descent.setEvaluationGoal(ev_goal);
    }
    void setGradientNormGoal(double norm_goal){
        gradient_descent.setGradientNormGoal(norm_goal);
    }
    SEXP getMinimalArgument(){
        return Rcpp::wrap<Purple::Vector<double>>(gradient_descent.getMinimalArgument());
    }
    void setInitialArgument(SEXP vect_inicial){
        gradient_descent.setInitialArgument(Rcpp::as<Purple::Vector<double>>(vect_inicial));
    }
private:
    Purple::SumaCuadrados funcion_objetivo;
    Purple::GradientDescent gradient_descent;
};

RCPP_MODULE(ModuloSumaCuadrados){
    using namespace Rcpp;

    class_<GradientDescentSumaCuadradosClase>( "GradientDescentSumaCuadradosClase" )

    .constructor<SEXP,SEXP>()

    .method( "solve", &GradientDescentSumaCuadradosClase::getMinimalArgument )
    .method( "setGradientNormGoal", &GradientDescentSumaCuadradosClase::setGradientNormGoal )
    .method( "setEvaluationGoal", &GradientDescentSumaCuadradosClase::setEvaluationGoal )
    .method( "starting_point", &GradientDescentSumaCuadradosClase::setInitialArgument )
    ;
}
