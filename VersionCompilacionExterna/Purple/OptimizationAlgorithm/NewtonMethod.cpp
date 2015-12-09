/******************************************************************************/
/*                                                                            */
/*   N E W T O N   M E T H O D   C L A S S                                    */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#include "NewtonMethod.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <time.h>

namespace Purple
{

// GENERAL CONSTRUCTOR 

/// General constructor. It creates a Newton's method object
/// associated to an objective function object.
/// It also initializes the class members to their default values:
///
/// Initial argument: Random point whithin the objective function domain.
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e69.
/// <li> Gradient norm goal: 0.0.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of iterations: 1000. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period: 25. 
/// </ul>
///
/// @param newObjectiveFunction: Pointer to an objective function object.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

NewtonMethod::NewtonMethod(ObjectiveFunction* newObjectiveFunction)
: OptimizationAlgorithm(newObjectiveFunction)
{
   // Optimization parameters

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> lowerBound = objectiveFunction->getLowerBound();
   Vector<double> upperBound = objectiveFunction->getUpperBound();

   Vector<double> newInitialArgument(numberOfVariables, 0.0);

   for(int i = 0; i < numberOfVariables; i++)
   {
      double random = (double)rand()/(RAND_MAX+1.0);

      newInitialArgument[i] 
      = lowerBound[i] + (upperBound[i] - lowerBound[i])*random;
   }

   initialArgument = newInitialArgument;

   // Stopping criteria

   evaluationGoal = 0.0;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 100;

   // User stuff

   showPeriod = 25;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a Newton's method optimization algorithm object
/// not associated to any objective function object.
/// It also initializes the class members to their default values:
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e69.
/// <li> Gradient norm goal: 0.0.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of iterations: 1000. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period: 25. 
/// </ul>
///
/// @see OptimizationAlgorithm.

NewtonMethod::NewtonMethod(void) : OptimizationAlgorithm()
{
   // Stopping criteria

   evaluationGoal = 0.0;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 100;

   // User stuff

   showPeriod = 25;
}


// DESTRUCTOR

/// Destructor.

NewtonMethod::~NewtonMethod(void)
{

}


// METHODS


// Vector<double> getInitialArgument(void)

/// This method returns the initial objective function argument to be used by 
/// the Newton's method for optimization. 

Vector<double> NewtonMethod::getInitialArgument(void)
{
   return(initialArgument);
}


// double getGradientNormGoal(void) method

/// This method returns the goal value for the norm of the objective function
/// gradient.
/// This is used as a stopping criterium when optimizing a function.

double NewtonMethod::getGradientNormGoal(void)
{
   return(gradientNormGoal);
}


// int getMaximumNumberOfIterations(void) method

/// This method returns the maximum number of iterations to be performed by the 
/// Newton's method during the optimization process. 
/// This is used as a stopping criterium when optimizing an objective function.

int NewtonMethod::getMaximumNumberOfIterations(void)
{
   return(maximumNumberOfIterations);
}


// int getShowPeriod(void) method

/// This method returns the number of iterations between the optimization 
/// showing progress. 

int NewtonMethod::getShowPeriod(void)
{
   return(showPeriod);    
}


// void setInitialArgument(Vector<double>) method

/// This method sets a new initial objective function argument to be used by 
/// the Newton's method for optimization. 
///
/// @param newInitialArgument: Initial argument Vector.

void NewtonMethod::setInitialArgument(Vector<double> newInitialArgument)
{
   int size = newInitialArgument.getSize();

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   if(size != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: NewtonMethod class. "
                << "double setInitialArgument(Vector<double>) method." << std::endl
                << "Size of initial argument must be equal to number of variables."
                << std::endl << std::endl;

      exit(1);
   }


   initialArgument = newInitialArgument;
}


// void setGradientNormGoal(double) method

/// This method sets a new the goal value for the norm of the 
/// objective function gradient. 
/// This is used as a stopping criterium when optimizing an objective function.
///
/// @param newGradientNormGoal: 
/// Goal value for the norm of the objective function gradient.

void NewtonMethod::setGradientNormGoal(double newGradientNormGoal)
{
   if(gradientNormGoal < 0.0)
   {
      std::cout << std::endl
                << "Error: NewtonMethod class." << std::endl
                << "void setGradientNormGoal(double) method."
                << std::endl
                << "Gradient norm goal must be equal or greater than 0."
                << std::endl << std::endl;
      exit(1);
   }

   // Set gradient norm goal

   gradientNormGoal = newGradientNormGoal;
}


// void setMaximumNumberOfIterations(int) method

/// This method sets a new maximum number of iterations in the optimization 
/// process. 
///
/// @param newMaximumNumberOfIterations: Maximum number of iterations.

void NewtonMethod
::setMaximumNumberOfIterations(int newMaximumNumberOfIterations)
{
   if(newMaximumNumberOfIterations <= 0)
   {
      std::cout << std::endl
                << "Error: NewtonMethod class." << std::endl
                << "void setMaximumNumberOfIterations(int) method."
                << std::endl
                << "Maximum number of iterations must be greater than 0."
                << std::endl
                << std::endl;

      exit(1);
   }

   // Set maximum number of iterations

   maximumNumberOfIterations = newMaximumNumberOfIterations;
}



// void setShowPeriod(int) method

/// This method sets a new number of iterations between the optimization
/// showing progress. 
///
/// @param newShowPeriod: Show period.

void NewtonMethod::setShowPeriod(int newShowPeriod)
{
   if(newShowPeriod <= 0)
   {
      std::cout << std::endl
                << "Error: NewtonMethod class." << std::endl
                << "void setShowPeriod(int) method."
                << std::endl
                << "Show period must be greater than 0."
                << std::endl << std::endl;

      exit(1);
   }

   // Set show period

   showPeriod = newShowPeriod;
}


// void getMinimalArgument(void) method

/// This method optimizes an objective function according to the 
/// Newton's method. 
/// It returns the minimal argument of the objective function.
/// Optimization occurs according to the optimization parameters. 

Vector<double> NewtonMethod::getMinimalArgument(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> minimalArgument(numberOfVariables, 0.0);
   Vector<double> argument(numberOfVariables, 0.0);

   // Evaluation history vector

   Vector<double> newEvaluationHistory(maximumNumberOfIterations+1, 0.0);

   evaluationHistory = newEvaluationHistory;

   // Gradient norm optimization history vector

   Vector<double> newGradientNormHistory(maximumNumberOfIterations+1, 0.0);

   gradientNormHistory = newGradientNormHistory;

   double evaluation = 0.0;

   Vector<double> gradient(numberOfVariables, 0.0);

   double gradientNorm = 0.0;

   Matrix<double> hessian(numberOfVariables, numberOfVariables, 0.0);
   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables, 0.0);
   Vector<double> inverseHessianGradientProduct(numberOfVariables, 0.0);

   time_t beginningTime, currentTime;

   double elapsedTime = 0.0;

   // Set beginning time 

   time(&beginningTime);

   // Main loop

   std::cout << std::endl
             << "Getting minimal argument with Newton's method..." 
             << std::endl;

   argument = initialArgument;
    
   // Initial objective function evaluation
   
   evaluation = objectiveFunction->getEvaluation(argument);

   evaluationHistory[0] = evaluation;

   if(evaluation <= evaluationGoal)
   {          
      std::cout << std::endl
                << "Initial evaluation is less than goal." << std::endl
                << "Initial evaluation: " << evaluation << std::endl;
      
      minimalArgument = argument;

      // Print minimal argument to screen

      std::cout << std::endl
                << "Minimal argument:" << std::endl;
   
      for(int i = 0; i < numberOfVariables; i++)
      {
         std::cout << minimalArgument[i] << " ";        
      }
      
      return(minimalArgument);        
   }
   else
   {
      std::cout << "Initial evaluation: " <<  evaluation << std::endl;      
   }

   // Initial objective function gradient

   gradient = objectiveFunction->getGradient(argument);

   gradientNorm = objectiveFunction->getGradientNorm(gradient);

   gradientNormHistory[0] = gradientNorm;

   if(gradientNorm <= gradientNormGoal)
   {          
      std::cout << std::endl
                << "Initial gradient norm is less than goal." << std::endl
                << "Initial gradient norm: " << gradientNorm << std::endl;
              
      minimalArgument = argument;
     
      // Print minimal argument to screen

      std::cout << std::endl
                << "Minimal argument:" << std::endl;
   
      for(int i = 0; i < numberOfVariables; i++)
      {
         std::cout << minimalArgument[i] << " ";        
      }
      
      return(minimalArgument);        
   }
   else
   {
      std::cout << "Initial gradient norm: " <<  gradientNorm << std::endl;      
   }

   // Loop over iterations

   for(int iteration = 1; iteration <= maximumNumberOfIterations; iteration++)
   {
      // Objective function Hessian

      inverseHessian = objectiveFunction->getInverseHessian(argument);

      // Inverse Hessian - gradient product

      for(int i = 0; i < numberOfVariables; i++)
      {
         for(int j = 0; j < numberOfVariables; j++)
         {
            inverseHessianGradientProduct[i] += inverseHessian[i][j]*gradient[j];
         }
      }

      // Get new argument

      for (int i = 0; i < numberOfVariables; i++)
      {
         argument[i] -= inverseHessianGradientProduct[i];
      }
      
      
      // Objective function evaluation
   
      evaluation = objectiveFunction->getEvaluation(argument);

      evaluationHistory[iteration] = evaluation;

      // Objective function gradient

      gradient = objectiveFunction->getGradient(argument);

      gradientNorm = objectiveFunction->getGradientNorm(gradient);

      gradientNormHistory[iteration] = gradientNorm;
      

      // Stopping Criteria

      // Evaluation goal 

      if (evaluation <= evaluationGoal)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Evaluation goal reached." << std::endl;

         std::cout << "Evaluation: " << evaluation  << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << std::endl;

         break;
      }

      // Norm of objective function gradient goal 

      if (gradientNorm <= gradientNormGoal)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Gradient norm goal reached."
                   << std::endl;  

         std::cout << "Evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;

         break;
      }

      // Maximum optimization time

      time(&currentTime);

      elapsedTime = difftime(currentTime, beginningTime);

      if (elapsedTime >= maximumTime)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Maximum optimization time reached."
                   << std::endl;

         std::cout << "Evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;

         break;
      }

      // Maximum number of iterations

      if (iteration == maximumNumberOfIterations)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Maximum number of iterations reached."
                   << std::endl;

         std::cout << "Evaluation: " << evaluation << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << std::endl;

         break;
      }

      // Progress

      if(iteration % showPeriod == 0)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << "; " << std::endl;

         std::cout << "Evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;
      }
   }

   // Set minimal argument

   minimalArgument = argument;

   // Print minimal argument to screen

   std::cout << std::endl
             << "Minimal argument:" << std::endl;
   
   for(int i = 0; i < numberOfVariables; i++)
   {
      std::cout << minimalArgument[i] << " ";        
   }

   std::cout << std::endl;
   
   return(minimalArgument);
}


// void print(void) method

/// This method prints to the screen the initial argumetn and the 
/// stopping criteria concerning the Newton's method object:
///
/// Initial argument.
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of iterations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period. 
/// </ul>

void NewtonMethod::print(void)
{
   std::cout << std::endl
             << "Newton's Method Object." << std::endl;

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   // Initial argument

   std::cout << std::endl
             << "Initial argument:" << std::endl;

   for(int i = 0; i < numberOfVariables; i++)
   {
      std::cout << initialArgument[i] << " ";        
   }
   
   std::cout << std::endl;

   // Stopping criteria

   std::cout << std::endl
             << "Stopping criteria: " << std::endl
             << "Evaluation goal: " << std::endl
             << evaluationGoal << std::endl
             << "Gradient norm goal" << std::endl 
             << gradientNormGoal <<std::endl
             << "Maximum time: " << std::endl
             << maximumTime << std::endl
             << "Maximum number of iterations: " << std::endl
             << maximumNumberOfIterations << std::endl;

   // User stuff

   std::cout << std::endl
             << "User stuff: " << std::endl
             << "Show period: " << std::endl
             << showPeriod
             << std::endl;

}


// void save(char*) method

/// This method saves the Newton's method object to a data file. 
///
/// Initial argument.
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of iterations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see load(char*).

void NewtonMethod::save(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: NewtonMethod class." << std::endl
                << "void save(char*) method."
                << std::endl
                << "Cannot open Newton method object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving Newton method object to data file..." << std::endl;
   }

   // Write file header

   file << "% Purple: An Open Source Numerical Optimization C++ Library." 
        << std::endl 
        << "% Newton Method Object." << std::endl; 

   int numberOfVariables = objectiveFunction->getNumberOfVariables();
 
   // Initial argument

   file << "InitialArgument:" << std::endl;

   for(int i = 0; i < numberOfVariables; i++)
   {
      file << initialArgument[i] << " ";        
   }
   
   file << std::endl;

   // Stopping criteria

   file << "EvaluationGoal:" << std::endl
        << evaluationGoal << std::endl
        << "GradientNormGoal:" << std::endl
        << gradientNormGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfIterations: " << std::endl
        << maximumNumberOfIterations << std::endl;

   // User stuff

   file << "ShowPeriod: " << std::endl
        << showPeriod << std::endl;

   file.close();
}


// void load(char*) method

/// This method loads a Newton method object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
///
/// Initial argument.
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Gradient norm goal.
/// <li> Maximum time.
/// <li> Maximum number of iterations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see save(char*).

void NewtonMethod::load(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Error: NewtonMethod class." << std::endl
                << "void load(char*) method."
                << std::endl
                << "Cannot open Newton method object data file."  << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading Newton method object from data file..."
                << std::endl;
   }

   std::string word;

   // Initial argument

   while(word != "InitialArgument:")
   {
      file >> word;
   }

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   for(int i = 0; i < numberOfVariables; i++)
   {
      file >> initialArgument[i];        
   }

   // Stopping criteria: 

   // Evaluation goal

   file >> word;

   file >> evaluationGoal;

   // Norm of objective function gradient goal

   file >> word;

   file >> gradientNormGoal;

   // Maximum time

   file >> word;

   file >> maximumTime;

   // Maximum number of iterations

   file >> word;

   file >> maximumNumberOfIterations;

   // User stuff: 

   // Iterations between showing progress

   file >> word;

   file >> showPeriod;

   // Close file

   file.close();
}


// void saveOptimizationHistory(char*) method 

/// This method saves the optimization history to a data file:
///
/// <ol>
/// <li> Iteration.
/// <li> Objective function evaluation.
/// <li> Objective function gradient norm.
/// </ol>
///
/// @param filename: Filename.

void NewtonMethod::saveOptimizationHistory(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   // Write file header 

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: NewtonMethod class. " << std::endl
                << "void saveOptimizationHistory(char*) method."
                << std::endl
                << "Cannot open optimization history data file." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl 
                << "Saving optimization history to data file..." << std::endl;
   }

   // Write file header

   file << "% Purple: An Open Source Numerical Optimization C++ Library." 
        << std::endl 
        << "% Newton Method Optimization History." << std::endl
        << "% 1 - Iteration." << std::endl
        << "% 2 - Objective function evaluation." << std::endl
        << "% 3 - Objective function gradient norm." << std::endl;

   // Write file data

   int size = evaluationHistory.getSize();

   for (int i = 0; i < size; i++)
   {
      file << i << " "
           << evaluationHistory[i] << " "
           << gradientNormHistory[i] << std::endl;
   }

   file << std::endl;

   file.close();
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
