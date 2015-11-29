/******************************************************************************/
/*                                                                            */
/*   N E W T O N   M E T H O D   C L A S S   H E A D E R                      */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#ifndef __NEWTONMETHOD_H__
#define __NEWTONMETHOD_H__

#include "OptimizationAlgorithm.h"
#include "../ObjectiveFunction/ObjectiveFunction.h"


namespace Purple
{

/// This concrete class represents the Newton's method optimization algorithm
/// for an objective function.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

class NewtonMethod : public OptimizationAlgorithm
{

private: 

   /// Initial argument

   Vector<double> initialArgument;

   /// Objective function gradient norm goal.
   /// It is used as a stopping criterion.

   double gradientNormGoal;

   /// Maximum number of iterations.
   /// It is used as a stopping criterion.

   int maximumNumberOfIterations;

   /// Number of iterations between the training showing progress.

   int showPeriod;

   /// Evaluation of objective function optimization history.

   Vector<double> evaluationHistory;

   /// Gradient norm of objective function optimization history.

   Vector<double> gradientNormHistory;


public:

   // GENERAL CONSTRUCTOR

   NewtonMethod(ObjectiveFunction*);


   // DEFAULT CONSTRUCTOR

   NewtonMethod(void);


   // DESTRUCTOR

   virtual ~NewtonMethod(void);


   // METHODS

   // Get methods

   Vector<double> getInitialArgument(void);

   double getGradientNormGoal(void);
   int getMaximumNumberOfIterations(void);

   int getShowPeriod(void);

   // Set methods

   void setInitialArgument(Vector<double>);

   void setGradientNormGoal(double);
   void setMaximumNumberOfIterations(int);

   void setShowPeriod(int);

   // Optimization methods

   Vector<double> getMinimalArgument(void);

   // Utility methods

   void print(void);

   void load(char*);
   void save(char*);

   void saveOptimizationHistory(char*);
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
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
