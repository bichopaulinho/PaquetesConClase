// Clase SumaCuadrados: función objetivo que se optimiza en un problema de regresión: Suma de cuadrados

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "SumaCuadrados.h"


namespace Purple
{

//constructores
SumaCuadrados::SumaCuadrados(void): ObjectiveFunction()
{

}

SumaCuadrados::SumaCuadrados(const Vector<double> y, const Matrix<double> X)
{
    if (X.getNumberOfRows() != y.getSize()){

        std::cout << std::endl
                  << "Error: SumaCuadrados class. "
                  << "Dimensiones y, X incoherentes" << std::endl
                  << std::endl;
        exit(1);
    }

    numberOfVariables = X.getNumberOfColumns();

    Matrix<double> ModelMatrix(X);
    Vector<double> TargetVariable(y);

}

//destructores
SumaCuadrados::~SumaCuadrados(void)
{

}

//getEvaluation

double SumaCuadrados::getEvaluation(Vector<double> argument)
{

    if (argument.getSize() != numberOfVariables){
        std::cout << std::endl
                  << "Error: SumaCuadrados class. Metodo getEvaluation " << std::endl
                  << "dimension vector coeficientes incorrecta" << std::endl
                  << std::endl;
        exit(1);
    }

    // valores estimados en la regresión
    Vector<double> PredictedValues(TargetVariable.getSize());

    // tbc

    // suma de cuadrados
    double ss=0.0;

    //tbc

    return ss;

}

// getModelMatrix
Matrix<double> SumaCuadrados::getModelMatrix(void)
{
    return(ModelMatrix);
}

// getTargetVariable
Vector<double> SumaCuadrados::getTargetVariable(void)
{
    return(TargetVariable);
}



}
