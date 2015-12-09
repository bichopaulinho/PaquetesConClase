/******************************************************************************/
/*                                                                            */
/*   O B J E C T I V E   F U N C T I O N   C L A S S   H E A D E R            */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __OBJECTIVEFUNCTION_H__
#define __OBJECTIVEFUNCTION_H__

#include "../Utilities/Vector.h"
#include "../Utilities/Matrix.h"

namespace Purple
{

/// This abstract class represents the concept of objective function.
/// Any derived class must implement the getEvaluation(Vector<double>) method.

class ObjectiveFunction
{

protected:

   /// Number of variables in the objective function.

   int numberOfVariables;

   /// Lower bound of objective function domain.

   Vector<double> lowerBound;

   /// Upper bound of objective function domain.

   Vector<double> upperBound;

   /// Number of calls to the getEvaluation(Vector<double>) method.
   ///
   /// @see getEvaluation(Vector<double>).

   int numberOfEvaluations;

   /// Epsilon value for numerical differentiation.

   double epsilon;

   // Utility methods
   
   double getDeterminant(Matrix<double>);

public:

   // GENERAL CONSTRUCTOR

   ObjectiveFunction(void);


   // DESTRUCTOR

   virtual ~ObjectiveFunction(void);


   // METHODS

   // Get methods

   int getNumberOfVariables(void);

   Vector<double> getLowerBound(void);
   Vector<double> getUpperBound(void);

   Matrix<double> getDomain(void);
   
   double getEpsilon(void);
   int getNumberOfEvaluations(void);

   // Set methods

   void setNumberOfVariables(int);

   void setLowerBound(Vector<double>);
   void setUpperBound(Vector<double>);

   void setDomain(Matrix<double>);

   void setEpsilon(double);
   void setNumberOfEvaluations(int);

   // Objective function methods

   /// This method returns the evaluation value of an objective function for
   /// a given argument.
   ///
   /// @see getGradient(Vector<double>).
   /// @see getHessian(Vector<double>).

   virtual double getEvaluation(Vector<double>) = 0;

   // Objective function gradient vector methods

   /// This method returns the objective function gradient vector for a
   /// given argument.
   ///
   /// @see getEvaluation(Vector<double>).
   /// @see getHessian(Vector<double>).

   virtual Vector<double> getGradient(Vector<double>);

   double getGradientNorm(Vector<double>);

   // Objective function Hessian matrix methods

   /// This method returns the objective function Hessian matrix for a
   /// given argument.
   ///
   /// @see getEvaluation(Vector<double>).
   /// @see getGradient(Vector<double>).

   virtual Matrix<double> getHessian(Vector<double>);

   /// This method returns the inverse of the Hessian matrix for a
   /// given argument.
   ///
   /// @see getEvaluation(Vector<double>).
   /// @see getGradient(Vector<double>).
   /// @see getHessian(Vector<double>).

   virtual Matrix<double> getInverseHessian(Vector<double>);

   // Utility methods

   /// This prints to the screen any useful information of the objective function
   /// for a given argument during the optimization process.

   virtual void print(Vector<double>);
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
