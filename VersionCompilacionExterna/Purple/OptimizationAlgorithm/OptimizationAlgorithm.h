/******************************************************************************/
/*                                                                            */
/*   O P T I M I Z A T I O N   A L G O R I T H M   C L A S S   H E A D E R    */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __OPTIMIZATIONALGORITHM_H__
#define __OPTIMIZATIONALGORITHM_H__

#include "../ObjectiveFunction/ObjectiveFunction.h"

namespace Purple
{

/// This abstract class represents the concept of optimization algorithm.
/// Any derived class must implement the getMinimalArgument(void) method
///
/// @see ObjectiveFunction.

class OptimizationAlgorithm
{

protected:

   // FIELDS

   /// Pointer to an objective function object.

   ObjectiveFunction* objectiveFunction;

   /// Objective function evaluation goal value. 
   /// It is used as a stopping criterion.

   double evaluationGoal;

   /// Maximum optimization time. It is used as a stopping criterion.

   double maximumTime;

public:


   // GENERAL CONSTRUCTOR

   OptimizationAlgorithm(ObjectiveFunction*);


   // DEFAULT CONSTRUCTOR

   OptimizationAlgorithm(void);


   // DESTRUCTOR

   virtual ~OptimizationAlgorithm(void);


   // METHODS

   // Get methods

   ObjectiveFunction* getObjectiveFunction(void);

   double getEvaluationGoal(void);
   double getMaximumTime(void);

   // Set methods

   void setObjectiveFunction(ObjectiveFunction*);
   
   void setEvaluationGoal(double);
   void setMaximumTime(double);

   // Optimization methods

   virtual Vector<double> getMinimalArgument(void) = 0;
};

}

#endif


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
