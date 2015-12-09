/******************************************************************************/
/*                                                                            */
/*   R A N D O M   S E A R C H   C L A S S   H E A D E R                      */
/*                                                                            */ 
/*   Roberto Lopez                                                            */ 
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */ 
/*                                                                            */
/******************************************************************************/

#ifndef __RANDOMSEARCH_H__
#define __RANDOMSEARCH_H__

#include "OptimizationAlgorithm.h"
#include "../ObjectiveFunction/ObjectiveFunction.h"
 
namespace Purple
{

/// This concrete class represents the random search optimization algorithm.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

class RandomSearch : public OptimizationAlgorithm
{

private: 

   // FIELDS

   /// Objective function evaluation history

   Vector<double> evaluationHistory;

   /// Maximum number of iterations.
   /// It is used as a train stopping criterion.

   int maximumNumberOfIterations;

   /// Number of iterations between the training showing progress.

   int showPeriod;

public:

   // GENERAL CONSTRUCTOR

   RandomSearch(ObjectiveFunction*); 


   // DEFAULT CONSTRUCTOR

   RandomSearch(void); 


   // DESTRUCTOR

   virtual ~RandomSearch(void);


   // METHODS

   // Get methods

   int getMaximumNumberOfIterations(void);
   int getShowPeriod(void);

   // Set methods

   void setMaximumNumberOfIterations(int);
   void setShowPeriod(int);

   // Optimization methods

   Vector<double> getMinimalArgument(void);

   // Utiltity methods

   void print(void);

   void load(char*);
   void save(char*);

   void saveOptimizationHistory(char*);
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
