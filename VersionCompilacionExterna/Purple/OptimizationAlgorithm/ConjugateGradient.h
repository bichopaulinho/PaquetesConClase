/******************************************************************************/
/*                                                                            */
/*   C O N J U G A T E   G R A D I E N T   C L A S S   H E A D E R            */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#ifndef __CONJUGATEGRADIENT_H__
#define __CONJUGATEGRADIENT_H__

#include "OptimizationAlgorithm.h"
#include "../ObjectiveFunction/ObjectiveFunction.h"

namespace Purple
{

/// This concrete class represents a conjugate gradient optimization algorithm
/// for an objective function.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

class ConjugateGradient : public OptimizationAlgorithm
{

public:

   // ENUMERATIONS

   /// Enumeration of the available optimization operators for the 
   /// search direction.

   enum SearchDirectionMethod{PolakRibiere, FletcherReeves};

   /// Enumeration of the available optimization operators for the optimal
   /// step sizes.

   enum OptimalStepSizeMethod{GoldenSection, BrentMethod};

private: 

   // FIELDS

   /// Initial argument

   Vector<double> initialArgument;

   /// Objective function gradient norm goal.
   /// It is used as a stopping criterion.

   double gradientNormGoal;

   /// Maximum number of iterations.
   /// It is used as an optimization stopping criterion.

   int maximumNumberOfIterations;

   /// Number of iterations between the optimization showing progress.

   int showPeriod;

   /// Initial step size in line minimization.

   double firstStepSize;

   /// Tolerance in optimal step size.

   double optimalStepSizeTolerance;

   /// Step size at wich a warning message is written to the screen during line
   /// minimization.

   double warningStepSize;

   /// Evaluation history.

   Vector<double> evaluationHistory;

   /// Norm of objective function gradient history.

   Vector<double> gradientNormHistory;

   // METHODS

   // Search direction methods

   double getPolakRibiereParameter(Vector<double>, Vector<double>);
   double getFletcherReevesParameter(Vector<double>, Vector<double>);

   Vector<double>
   getFletcherReevesSearchDirection(Vector<double>, Vector<double>, Vector<double>);

   Vector<double>
   getPolakRibiereSearchDirection(Vector<double>, Vector<double>, Vector<double>);

   // Optimal step size methods

   double getGoldenSectionOptimalStepSize
   (double, double, Vector<double>, Vector<double>);
   
   double getBrentMethodOptimalStepSize
   (double, double, Vector<double>, Vector<double>);

   // Search direction optimization operators enumeration.

   SearchDirectionMethod searchDirectionMethod;

   // Optimal step size optimization operators enumeration.

   OptimalStepSizeMethod optimalStepSizeMethod;

   // Utility methods

   double getMinimum(Vector<double>);
   double getMaximum(Vector<double>);


public:

   // GENERAL CONSTRUCTOR

   ConjugateGradient(ObjectiveFunction*);


   // DEFAULT CONSTRUCTOR

   ConjugateGradient(void); 


   // DESTRUCTOR

   virtual ~ConjugateGradient(void);


   // METHODS

   // Get methods

   SearchDirectionMethod getSearchDirectionMethod(void);
   OptimalStepSizeMethod getOptimalStepSizeMethod(void);
   Vector<double> getInitialArgument(void);
   double getGradientNormGoal(void);
   int getMaximumNumberOfIterations(void);
   int getShowPeriod(void);
   double getFirstStepSize(void);
   double getOptimalStepSizeTolerance(void);
   double getWarningStepSize(void);

   // Set methods

   void setSearchDirectionMethod(SearchDirectionMethod);
   void setOptimalStepSizeMethod(OptimalStepSizeMethod);
   void setFirstStepSize(double);
   void setOptimalStepSizeTolerance(double);
   void setInitialArgument(Vector<double>);
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
