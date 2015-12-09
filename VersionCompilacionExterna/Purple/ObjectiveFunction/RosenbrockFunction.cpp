/******************************************************************************/
/*                                                                            */
/*   R O S E N B R O C K   F U N C T I O N   C L A S S                        */
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

#include "RosenbrockFunction.h"

namespace Purple
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a Rosenbrock's objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables = 2.
/// <li> Lower bound = -5.12,...,-5.12.
/// <li> Upper bound = 5.12,...,5.12.
/// </ul> 

RosenbrockFunction::RosenbrockFunction(void) : ObjectiveFunction()
{
   numberOfVariables = 2;

   Vector<double> newLowerBound(numberOfVariables, -1.0);

   lowerBound = newLowerBound;

   Vector<double> newUpperBound(numberOfVariables, 1.0);
   
   upperBound = newUpperBound;
}


// DESTRUCTOR

/// Destructor.

RosenbrockFunction::~RosenbrockFunction(void)
{

}


// METHODS

// double getEvaluation(Vector<double>) method

/// This method returns the Rosenbrock's function evaluation for a given argument.
///
/// @param argument: Objective function argument.

double RosenbrockFunction::getEvaluation(Vector<double> argument)
{
   double evaluation = 0.0;

   int size = argument.getSize();

   if(size != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: RosenbrockFunction class. "
                << "double getEvaluation(Vector<double>) method." << std::endl
                << "Size of argument must be equal to number of variables." << std::endl
                << std::endl;

      std::cout << "Size of argument: " << size << std::endl
                << "Number of variables: " << numberOfVariables <<std::endl;

      exit(1);
   }

   // Get evaluation

   for(int i = 0; i < numberOfVariables-1; i++)
   {
      evaluation += 100.0*pow(argument[i+1] - pow(argument[i],2), 2)
      + pow(1.0 - argument[i], 2);
   }

   return(evaluation);
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
