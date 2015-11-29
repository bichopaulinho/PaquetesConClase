/******************************************************************************/
/*                                                                            */
/*   D E   J O N G   F U N C T I O N   C L A S S                              */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#include <iostream>
#include <fstream>
#include <math.h>

#include "DeJongFunction.h"

namespace Purple
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a De Jong's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables = 2.
/// <li> Lower bound = -5.12,...,-5.12.
/// <li> Upper bound = 5.12,...,5.12.
/// </ul> 

DeJongFunction::DeJongFunction(void) : ObjectiveFunction()
{
   numberOfVariables = 2;
   
   Vector<double> newLowerBound(numberOfVariables, -5.12);

   lowerBound = newLowerBound;

   Vector<double> newUpperBound(numberOfVariables, 5.12);
   
   upperBound = newUpperBound;
}


// DESTRUCTOR

/// Destructor.

DeJongFunction::~DeJongFunction(void)
{

}


// METHODS

// double getEvaluation(Vector<double>) method

/// This method returns the De Jong's function evaluation for a given argument.
///
/// @param argument: Objective function argument.

double DeJongFunction::getEvaluation(Vector<double> argument)
{
   double evaluation = 0.0;

   int size = argument.getSize();

   if(size != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: DeJongFunction class. "
                << "double getEvaluation(Vector<double>) method." << std::endl
                << "Size of argument must be equal to number of variables." << std::endl
                << std::endl;

      exit(1);
   }

   for(int i = 0; i < numberOfVariables; i++)
   {
      evaluation += pow(argument[i], 2);
   }

   return(evaluation);
}


// Vector<double> getGradient(Vector<double>) method

/// This method returns the De Jong's analytical gradient Vector 
/// for a given argument.
///
/// @param argument: Point at which the gradient is to be computed.

Vector<double> DeJongFunction::getGradient(Vector<double> argument)
{
   Vector<double> gradient(numberOfVariables, 0.0);

   for(int i = 0; i < numberOfVariables; i++)
   {
      gradient[i] = 2.0*argument[i];
   }

   return(gradient);
}


// Matrix<double> getHessian(Vector<double>) method

/// This method returns the De Jong's analytical Hessian Matrix 
/// for a given argument.
///
/// @param argument: Point at which the Hessian is to be computed.

Matrix<double> DeJongFunction::getHessian(Vector<double> argument)
{
   Matrix<double> hessian(numberOfVariables, numberOfVariables, 0.0);

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            hessian[i][j] = 2.0;
         }
         else
         {
            hessian[i][j] = 0.0;
         }
      }
   }

   return(hessian);
}


// Matrix<double> getInverseHessian(Vector<double>) method

/// This method returns the De Jong's analytical inverse Hessian Matrix 
/// for a given argument.
///
/// @param argument: Point at which the inverse Hessian is to be computed.

Matrix<double> DeJongFunction::getInverseHessian(Vector<double> argument)
{
   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables, 0.0);

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         if(i == j)
         {
            inverseHessian[i][j] = 0.5;
         }
         else
         {
            inverseHessian[i][j] = 0.0;
         }
      }
   }

   return(inverseHessian);
}

}


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
