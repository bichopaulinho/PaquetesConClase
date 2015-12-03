
#ifndef __SUMA_CUADRADOS__
#define __SUMA_CUADRADOS__

#include "ObjectiveFunction.h"


namespace Purple
{

/// Clase que representa la suma de cuadrados de una regresión
///
/// @see ObjectiveFunction.

class SumaCuadrados : public ObjectiveFunction
{

protected:
    Matrix<double> ModelMatrix;
    Vector<double> TargetVariable;

public:

    // GENERAL CONSTRUCTOR

    SumaCuadrados(void);

    // Constructor a partir del vector de la variable predicha y la matriz del modelo
    SumaCuadrados(Vector<double>, Matrix<double>);

    // DESTRUCTOR

    virtual ~SumaCuadrados(void);


    // METHODS

    // Objective function methods

    double getEvaluation(Vector<double>);

    Matrix<double> getModelMatrix(void);
    Vector<double> getTargetVariable(void);


};
}
#endif
