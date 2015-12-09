// Clase SumaCuadrados: función objetivo que se optimiza en un problema de regresión: Suma de cuadrados

#include <iostream>
#include <fstream>
#include <math.h>
#include <stdlib.h>

#include "ObjectiveFunction/SumaCuadrados.h"


namespace Purple
{

//constructores
SumaCuadrados::SumaCuadrados(void): ObjectiveFunction()
{

}

SumaCuadrados::SumaCuadrados(const Vector<double> y, const Matrix<double> X) : ModelMatrix(X),TargetVariable(y)
{
    if (X.getNumberOfRows() != y.getSize()){

        std::cout << std::endl
                  << "Error: SumaCuadrados class. "
                  << "Dimensiones y, X incoherentes" << std::endl
                  << std::endl;
        exit(1);
    }

    numberOfVariables = X.getNumberOfColumns();

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

    Vector<double> res = Residuals(argument);
    // tbc


    // error = sum(error^2)
    double error = 0;
    for (int irow=0;irow<ModelMatrix.getNumberOfRows();irow++){
    	error+=res[irow]*res[irow];
    }

    //tbc
    return error;
}
Vector<double> SumaCuadrados::Residuals(const Vector<double> & argument){
    // error = (ModelMatrix * vector - TargetVariable);
	Vector<double> res(TargetVariable.getSize());

    for(int irow = 0;irow<ModelMatrix.getNumberOfRows();irow++){
    	res[irow] = -TargetVariable[irow];
    	for(int jcol=0;jcol<ModelMatrix.getNumberOfColumns();jcol++){
    		res[irow]+=ModelMatrix[irow][jcol]*argument[jcol];
    	}
    }
    return res;

}

// getModelMatrix
Matrix<double> SumaCuadrados::getModelMatrix(void)
{
    return(ModelMatrix);
}

Vector<double> SumaCuadrados::getGradient(Vector<double> argument) {
	//return 2*transpose(ModelMatrix)*(res);
	Vector<double> salida(argument.getSize());
	Vector<double> res = Residuals(argument);

	for(int icol = 0;icol<ModelMatrix.getNumberOfColumns();icol++){
    	salida[icol] = 0;
    	for(int jrow=0;jrow<ModelMatrix.getNumberOfRows();jrow++){
    		salida[icol]+=2*ModelMatrix[jrow][icol]*res[jrow];
    	}
    }
	return salida;
}

// getTargetVariable
Vector<double> SumaCuadrados::getTargetVariable(void)
{
    return(TargetVariable);
}



}
