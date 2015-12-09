/******************************************************************************/
/*                                                                            */
/*   R A S T R I G I N   F U N C T I O N   C L A S S   H E A D E R            */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __RASTRIGINFUNCTION_H__
#define __RASTRIGINFUNCTION_H__

#include "ObjectiveFunction.h"

namespace Purple
{

/// This class represents the Rastrigin's objective function.
///
/// @see ObjectiveFunction.

class RastriginFunction : public ObjectiveFunction
{

public:

   // GENERAL CONSTRUCTOR

   RastriginFunction(void);


   // DESTRUCTOR

   virtual ~RastriginFunction(void);


   // METHODS

   // Objective function methods

   double getEvaluation(Vector<double>);

   // Objective function gradient vector methods

   Vector<double> getGradient(Vector<double>);

   // Objective function Hessian matrix methods

   Matrix<double> getHessian(Vector<double>);
   
   Matrix<double> getInverseHessian(Vector<double>);
};

}

#endif


// Purple: An Open Source Numerical Optimization C++ Library.
// Copyright (C) 2006 Roberto Lopez
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
