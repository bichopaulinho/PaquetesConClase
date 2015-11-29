/******************************************************************************/
/*                                                                            */
/*   O P T I M I Z A T I O N   A L G O R I T H M   C L A S S                  */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#include<iostream>

#include "OptimizationAlgorithm.h"

namespace Purple
{

// GENERAL CONSTRUCTOR
//
/// General constructor. It creates an optimization algorithm object associated
/// to an objective function object. 
///
/// @param newObjectiveFunction: Pointer to an objective function object.
///
/// @see ObjectiveFunction.

OptimizationAlgorithm
::OptimizationAlgorithm(ObjectiveFunction* newObjectiveFunction)
{
   objectiveFunction = newObjectiveFunction;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates an optimization algorithm object not associated
/// to any objective function object.

OptimizationAlgorithm::OptimizationAlgorithm(void)
{
   objectiveFunction = NULL;
}


// DESTRUCTOR 

/// Destructor.

OptimizationAlgorithm::~OptimizationAlgorithm(void)
{ 

}


// METHODS

// ObjectiveFunction* getObjectiveFunction(void) method

/// This method returns a pointer to the objective function object to which
/// the optimization algorithm is associated.

ObjectiveFunction* OptimizationAlgorithm::getObjectiveFunction(void)
{
   return(objectiveFunction);
}


// double getEvaluationGoal(void) method

/// This method returns the objective function evaluation goal value.
/// This is used as a stopping criterium when optimizing a function.
///  
/// @see getMinimalArgument(void).

double OptimizationAlgorithm::getEvaluationGoal(void)
{
   return(evaluationGoal);
}


// int getMaximumTime(void) method

/// This method returns the maximum optimization time.
///  
/// @see getMinimalArgument(void).

double OptimizationAlgorithm::getMaximumTime(void)
{
   return(maximumTime);
}



// void setObjectiveFunction(ObjectiveFunction*) method

/// This method sets a pointer to an objective function object to be associated 
/// to the optimization algorithm.
///
/// @param newObjectiveFunction: Pointer to an objective function object. 
/// 
/// @see getMinimalArgument(void).

void OptimizationAlgorithm
::setObjectiveFunction(ObjectiveFunction* newObjectiveFunction)
{
   objectiveFunction = newObjectiveFunction;
}


// void setEvaluationGoal(double) method

/// This method sets a new goal value for the objective function evaluation. 
/// This is used as a stopping criterium when optimizing an objective function.
///
/// @param newEvaluationGoal: Goal value for the evaluation.
/// 
/// @see getMinimalArgument(void).

void OptimizationAlgorithm::setEvaluationGoal(double newEvaluationGoal)
{
   evaluationGoal = newEvaluationGoal;
}


// void setMaximumTime(double) method

/// This method sets a new maximum optimization time.  
///
/// @param newMaximumTime: Maximum optimization time.
/// 
/// @see getMinimalArgument(void).

void OptimizationAlgorithm::setMaximumTime(double newMaximumTime)
{
   if(newMaximumTime <= 0.0)
   {
      std::cout << std::endl 
                << "Error: OptimizationAlgorithm class. " << std::endl
                << "void setMaximumTime(double) method." << std::endl
                << "Maximum time must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set maximum time
   
   maximumTime = newMaximumTime;
}

}


// Purple: An Open Source Numerical Optimization C++ Library.
// Copyright (C) 2005 Roberto Lopez 
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
