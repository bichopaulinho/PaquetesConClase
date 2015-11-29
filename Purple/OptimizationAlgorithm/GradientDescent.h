/******************************************************************************/
/*                                                                            */
/*   G R A D I E N T   D E S C E N T   C L A S S   H E A D E R                */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#ifndef __GRADIENTDESCENT_H__
#define __GRADIENTDESCENT_H__

#include "OptimizationAlgorithm.h"
#include "../ObjectiveFunction/ObjectiveFunction.h"


namespace Purple
{

/// This concrete class represents the gradient descent optimization algorithm
/// for an objective function.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

class GradientDescent : public OptimizationAlgorithm
{

public:

   // ENUMERATIONS

   /// Available optimization operators for obtaining the optimal step size.

   enum OptimalStepSizeMethod{GoldenSection, BrentMethod};


private: 

   /// Initial argument

   Vector<double> initialArgument;

   /// Objective function gradient norm goal.
   /// It is used as a stopping criterion.

   double gradientNormGoal;

   /// Maximum number of iterations.
   /// It is used as a stopping criterion.

   int maximumNumberOfIterations;

   /// Number of iterations between the optimization showing progress.

   int showPeriod;

   /// Inititial step size in line search for first iteration of gradient descent.

   double firstStepSize;

   /// Tolerance for the optimal step size.

   double optimalStepSizeTolerance;

   /// Step size value at wich a warning message is written to the screen.

   double warningStepSize;

   /// Evaluation of objective function optimization history.

   Vector<double> evaluationHistory;

   /// Gradient norm of objective function optimization history.

   Vector<double> gradientNormHistory;

   // METHODS

   // Optimal step size methods

   double getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>);
   double getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>);

   // Optimal step size methods enumeration

   OptimalStepSizeMethod optimalStepSizeMethod;

   // Utility methods

   double getMinimum(Vector<double>);
   double getMaximum(Vector<double>);


public:

   // GENERAL CONSTRUCTOR

   GradientDescent(ObjectiveFunction*);


   // DEFAULT CONSTRUCTOR

   GradientDescent(void); 


   // DESTRUCTOR

   virtual ~GradientDescent(void);


   // METHODS

   // Get methods

   Vector<double> getInitialArgument(void);

   OptimalStepSizeMethod getOptimalStepSizeMethod(void);

   double getFirstStepSize(void);
   double getOptimalStepSizeTolerance(void);
   
   double getGradientNormGoal(void);
   int getMaximumNumberOfIterations(void);

   int getShowPeriod(void);
   double getWarningStepSize(void);

   // Set methods

   void setInitialArgument(Vector<double>);

   void setOptimalStepSizeMethod(OptimalStepSizeMethod);

   void setFirstStepSize(double);
   void setOptimalStepSizeTolerance(double);

   void setGradientNormGoal(double);
   void setMaximumNumberOfIterations(int);

   void setShowPeriod(int);
   void setWarningStepSize(double);

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

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
