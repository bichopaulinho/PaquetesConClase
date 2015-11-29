/******************************************************************************/
/*                                                                            */
/*   P L A N E - C Y L I N D E R   C L A S S                                  */
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

#include "PlaneCylinder.h"

namespace Purple
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a plane-cylinder objective function object.
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of variables = 2.
/// <li> Lower bound = -5.12,...,-5.12.
/// <li> Upper bound = 5.12,...,5.12.
/// <li> Penalty = 100.
/// </ul> 

PlaneCylinder::PlaneCylinder(void) : ObjectiveFunction()
{
   numberOfVariables = 2;

   Vector<double> newLowerBound(numberOfVariables, -5.12);

   lowerBound = newLowerBound;

   Vector<double> newUpperBound(numberOfVariables, 5.12);
   
   upperBound = newUpperBound;

   penalty = 100.0;
}


// DESTRUCTOR

/// Destructor.

PlaneCylinder::~PlaneCylinder(void)
{

}


// METHODS

// double getPenalty(void) method

/// This method returns the penalty term ratio to be used in the plane-cylinder
/// problem
///
/// @see getError(Vector<double>)
/// @see getEvaluation(Vector<double>)

double PlaneCylinder::getPenalty(void)
{
   return(penalty);
}


/// This method sets a new penalty term ratio to be used in the plane-cylinder
/// problem
///
/// @param newPenalty: New penalty term ratio.
///
/// @see getError(Vector<double>)
/// @see getEvaluation(Vector<double>)

void PlaneCylinder::setPenalty(double newPenalty)
{
   penalty = newPenalty;
}


// double getError(Vector<double>) method

/// This method returns the error made in the constraint by a given argument. 
///
/// @param argument: argument.
///
/// @see getEvaluation(Vector<double>)

double PlaneCylinder::getError(Vector<double> argument)
{
   double error = 0.0;

   double x = argument[0];
   double y = argument[1];
   
   if(pow(x,2) + pow(y,2) <= 1.0)
   {
      error = 0.0;
   }
   else
   {
      error = pow(x,2) + pow(y,2) - 1.0;
   } 
   
   return(error);       
}


// double getEvaluation(void) method

/// This method returns the plane-cylinder function evaluation for a given 
/// argument.
///
/// @param argument: Objective function argument.

double PlaneCylinder::getEvaluation(Vector<double> argument)
{
   double evaluation = 0.0;

   int size = argument.getSize();

   if(size != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: PlaneCylinder class. "
                << "double getEvaluation(Vector<double>) method." << std::endl
                << "Size of argument must be equal to number of variables." << std::endl
                << std::endl;

      exit(1);
   }

   double x = argument[0];
   double y = argument[1];

   double error = getError(argument);

   evaluation = x + y + penalty*pow(error, 2);

   return(evaluation);
}


// void print(void) method

/// This method prints to the screen the error made in the constraint by 
/// a given argument during the optimization process.  
///
/// @param argument: argument.
///
/// @see getError(Vector<double>)

void PlaneCylinder::print(Vector<double> argument)
{
   double error = getError(argument);
   
   std::cout << "Error: " << error << std::endl;
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
