/******************************************************************************/
/*                                                                            */
/*   R A N D O M   S E A R C H   C L A S S                                    */
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
#include <algorithm>
#include <functional>
#include <math.h>
#include <time.h>

#include "RandomSearch.h"


namespace Purple
{

// GENERAL CONSTRUCTOR 

/// General constructor. It creates a random search optimization
/// algorithm object associated to an objective function object. 
/// It also initializes the class members to their default values:
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e69.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of iterations: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period: 10. 
/// </ul>
///
/// @param newObjectiveFunction: Pointer to an objective function object.

RandomSearch::RandomSearch(ObjectiveFunction* newObjectiveFunction)
: OptimizationAlgorithm(newObjectiveFunction)
{   
   // Stopping criteria

   evaluationGoal = -1.0e69;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 100; 
   
   // User stuff
   
   showPeriod = 10;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a random search optimization algorithm 
///object not associated to any objective function object. 
/// It also initializes the class members to their default values:
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e69.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of evaluations: 100. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Show period: 10. 
/// </ul>

RandomSearch::RandomSearch(void) : OptimizationAlgorithm()
{
   // Stopping criteria

   evaluationGoal = -1.0e69;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 100; 
   
   // User stuff
   
   showPeriod = 10;
}


// DESTRUCTOR

/// Destructor.

RandomSearch::~RandomSearch(void)
{

}

// METHODS

// int getMaximumNumberOfIterations(void) method

/// This method returns the maximum number of iterations in the optimization
/// process.

int RandomSearch::getMaximumNumberOfIterations(void)
{
   return(maximumNumberOfIterations);
}


// int getShowPeriod(void) method

/// This method returns the number of iterations between the optimization 
/// showing progress. 

int RandomSearch::getShowPeriod(void)
{
   return(showPeriod);    
}


// void setMaximumNumberOfIterations(int) method

/// This method sets a maximum number of iterations for optimization. 
/// Each iteration with the random search optimization algorithms costs one 
/// objective function evaluation.
///
/// @param newMaximumNumberOfIterations: 
/// Maximum number of iterations for optimization.

void RandomSearch
::setMaximumNumberOfIterations(int newMaximumNumberOfIterations)
{
   if(newMaximumNumberOfIterations <= 0)
   {
      std::cout << std::endl
                << "Error: RandomSearch Class." << std::endl
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

void RandomSearch::setShowPeriod(int newShowPeriod)
{
   if(newShowPeriod <= 0)
   {
      std::cout << std::endl
                << "Error: RandomSearch class." << std::endl
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

/// This method returns the miminal argument of the associated 
/// objective function according to the random search optimization algorithm.
/// Optimization occurs according to the optimization parameters. 

Vector<double> RandomSearch::getMinimalArgument(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();	

   Vector<double> minimalArgument(numberOfVariables, 0.0);

   double bestEvaluation = 0.0;
   
   Vector<double> argument(numberOfVariables, 0.0);
        
   double evaluation = 0.0;
   
   // Evaluation history vector

   Vector<double> newEvaluationHistory(maximumNumberOfIterations, 0.0);

   evaluationHistory = newEvaluationHistory;

   // Objective function domain   
   
   Vector<double> lowerBound = objectiveFunction->getLowerBound();
   Vector<double> upperBound = objectiveFunction->getUpperBound();
    
   time_t beginningTime, currentTime;

   double elapsedTime = 0.0;

   // Set beginning time 

   time(&beginningTime);

   std::cout << std::endl
             << "Getting minimal argument with random search..." 
             << std::endl;
 
   // Initial objective function evaluation

   for(int i = 0; i < numberOfVariables; i++)
   {
      double random = (double)rand()/(RAND_MAX+1.0);

      argument[i] = lowerBound[i] + (upperBound[i]- lowerBound[i])*random;
   }
      
   evaluation = objectiveFunction->getEvaluation(argument);

   // Set best evaluation and minimal argument

   bestEvaluation = evaluation;
   minimalArgument = argument;

   evaluationHistory[0] = evaluation;

   if(bestEvaluation <= evaluationGoal)
   {          
      std::cout << std::endl
                << "Initial evaluation is less than goal." << std::endl
                << "Initial evaluation: " << bestEvaluation << std::endl;

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
   else
   {
      std::cout << "Initial evaluation: " <<  bestEvaluation << std::endl;      
   }

   // Main loop

   for(int iteration = 1; iteration < maximumNumberOfIterations; iteration++)
   {
      for(int i = 0; i < numberOfVariables; i++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         argument[i] = lowerBound[i] + (upperBound[i]- lowerBound[i])*random;
      }

      evaluation = objectiveFunction->getEvaluation(argument);

      if(evaluation < bestEvaluation)
      {
         // Set best evaluation              
                    
         bestEvaluation = evaluation;

         // Set minimal argument

         minimalArgument = argument;         
      }

      evaluationHistory[iteration] = bestEvaluation;

      // Stopping Criteria

      // Evaluation goal 

      if (evaluation <= evaluationGoal)
      {
         std::cout << "Iteration " << iteration << ": " 
                   << "Evaluation goal reached." << std::endl; 

         std::cout << "Final evaluation: " << bestEvaluation << ";" << std::endl;

         objectiveFunction->print(argument);

         break;
      }

      // Maximum time

      time(&currentTime);

      elapsedTime = difftime(currentTime, beginningTime);

      if (elapsedTime >= maximumTime)
      {
         std::cout << "Iteration " << iteration << ": "
                   << "Maximum time reached." 
                   << std::endl;

         std::cout << "Final evaluation: " << bestEvaluation << ";" << std::endl;

         objectiveFunction->print(argument);

         break;
      }

      // Maximum number of iterations

      if (iteration == maximumNumberOfIterations-1)
      {
         std::cout << "Iteration " << iteration << ": " 
                   << "Maximum number of iterations reached." 
                   << std::endl;

         std::cout << "Final evaluation: " << bestEvaluation << ";" << std::endl;

         objectiveFunction->print(argument);

         break;
      }

      // Progress

      if(iteration % showPeriod == 0) 			
      {
         std::cout << "Iteration " << iteration << "; " << std::endl
                   << "Best evaluation: " << bestEvaluation << std::endl;

         objectiveFunction->print(argument);
      }
   }

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

/// This method prints to the screen 
/// the optimization parameters concerning the random search object:
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

void RandomSearch::print(void)
{
   std::cout << std::endl
             << "RandomSearch Object." << std::endl;

   // Stopping criteria

   std::cout << std::endl
             << "Stopping criteria: " << std::endl
             << "Evaluation goal: " << std::endl
             << evaluationGoal << std::endl
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

/// This method saves the random search object to a data file. 
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of evaluations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Evaluations between showing progress. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see load(char*).

void RandomSearch::save(char* filename)
{

   int numberOfVariables = objectiveFunction->getNumberOfVariables();	

   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: RandomSearch class." << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open random search object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving random search object to data file..." << std::endl;
   }

   // Write file header

   file << "% Purple: An Open Source Numerical Optimization C++ Library." 
        << std::endl 
        << "% Random Search Object." << std::endl; 

   // Stopping criteria

   file << "EvaluationGoal:" << std::endl
        << evaluationGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfEvaluations: " << std::endl
        << maximumNumberOfIterations << std::endl;

   // User stuff

   file << "ShowPeriod: " << std::endl
        << showPeriod << std::endl;

   file.close();
}


// void load(char*) method

/// This method loads a random search object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of evaluations. 
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Evaluations between showing progress. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see save(char*).

void RandomSearch::load(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Error: RandomSearch class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open random search object data file."  << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading random seach object from data file..."  << std::endl;
   }

   std::string word;


   // Stopping criteria: 

   // Initial argument

   while(word != "EvaluationGoal:")
   {
      file >> word;
   }

   file >> evaluationGoal;

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

/// This method saves the optimization history to a data file. 
///
/// <ol>
/// <li> Iteration.
/// <li> Objective function evaluation.
/// </ol>
///
/// @param filename: Filename.

void RandomSearch::saveOptimizationHistory(char* filename)
{
   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: RandomSearch class. "
                << "void saveOptimizationHistory(char*) method." << std::endl
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
        << "% Gradient Descent Optimization History." << std::endl
        << "% 1 - Iteration." << std::endl
        << "% 2 - Objective function evaluation." << std::endl;

   // Write file data

   int size = evaluationHistory.getSize();

   for (int i = 0; i < size; i++)
   {
      file << i << " "
           << evaluationHistory[i] << std::endl;
   }

   file << std::endl;

   file.close();
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
