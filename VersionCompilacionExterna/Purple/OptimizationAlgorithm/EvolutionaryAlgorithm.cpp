/******************************************************************************/
/*                                                                            */
/*   E V O L U T I O N A R Y   A L G O R I T H M   C L A S S                  */
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

#include "EvolutionaryAlgorithm.h"

namespace Purple
{

// GENERAL CONSTRUCTOR 
//
/// General constructor. It creates an evolutionary optimization
/// algorithm object associated to an objective function object. 
/// Population is initialized at random within the objective function domain.
/// The constructor also initializes the class members to their default values:
///
/// Optimization operators:
/// <ul>
/// <li> Fitness assignment method: Linear ranking.
/// <li> Selection method: Stochastic universal sampling.
/// <li> Recombination method: Intermediate.
/// <li> Mutation method: Normal.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> Population size: 10 * [number of free parameters].
/// <li> Selective pressure: 1.5.
/// <li> Recombination size: 0.25.
/// <li> Mutation rate: = 1.0 / [number of free parameters].
/// <li> Mutation range: = 0.1
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: 0.0.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of generations: 100.
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Generations between showing progress: 10.
/// </ul>
///
/// @param newObjectiveFunction: Pointer to an objective function object.

EvolutionaryAlgorithm::EvolutionaryAlgorithm(ObjectiveFunction* newObjectiveFunction)
: OptimizationAlgorithm(newObjectiveFunction)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   // Fitness assignment method

   fitnessAssignmentMethod = LinearRanking;

   // Selection method

   selectionMethod = RouletteWheel;

   // Recombination method

   recombinationMethod = Intermediate;

   // Mutation method

   mutationMethod = Normal;

   // Optimization parameters

   if(numberOfVariables == 0)
   {                     
      populationSize = 0;

      selectivePressure = 0.0;

      recombinationSize = 0.0;

      mutationRate = 0.0;
      mutationRange = 0.0;
   }
   else
   {
      populationSize = 10*numberOfVariables;

      selectivePressure = 1.5;

      recombinationSize = 0.25;

      mutationRate = 1.0/numberOfVariables;
      mutationRange = 0.1;       
   }

   // Stopping criteria

   evaluationGoal = -1.0e69;
   maximumTime = 1.0e69;

   maximumNumberOfGenerations = 100;

   showPeriod = 1;

   // Population matrix

   Vector<double> lowerBound = objectiveFunction->getLowerBound();
   Vector<double> upperBound = objectiveFunction->getUpperBound();

   Matrix<double> newPopulation(populationSize, numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         newPopulation[i][j]  
         = lowerBound[j] + (upperBound[j] - lowerBound[j])*random;         
      }
   }

   population = newPopulation;

   // Evaluation vector

   Vector<double> newEvaluation(populationSize, 0.0);

   evaluation = newEvaluation;

   // Fitness vector

   Vector<double> newFitness(populationSize, 0.0);

   fitness = newFitness;

   // Selection vector

   Vector<bool> newSelection(populationSize, false);

   selection = newSelection;
}


// DEFAULT CONSTRUCTOR
//
/// Default constructor. It creates an evolutionary algorithm object
/// not associated to any objective function object. 
/// It also initializes the class members to their default values:
///
/// Optimization operators:
/// <ul>
/// <li> Fitness assignment method: Linear ranking.
/// <li> Selection method: Stochastic universal sampling.
/// <li> Recombination method: Intermediate.
/// <li> Mutation method: Normal.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> Population size: 0.
/// <li> Selective pressure: 0.0.
/// <li> Recombination size: 0.0.
/// <li> Mutation rate: = 0.0.
/// <li> Mutation range: = 0.0.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal: -1.0e69.
/// <li> Maximum time: 1.0e6.
/// <li> Maximum number of generations: 100.
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Generations between showing progress: 10.
/// </ul>

EvolutionaryAlgorithm::EvolutionaryAlgorithm(void) : OptimizationAlgorithm()
{
   // Fitness assignment method

   fitnessAssignmentMethod = LinearRanking;

   // Selection method

   selectionMethod = StochasticUniversalSampling;

   // Recombination method

   recombinationMethod = Intermediate;

   // Mutation method

   mutationMethod = Normal;

   // Optimization parameters

   populationSize = 0;

   selectivePressure = 0.0;

   recombinationSize = 0.0;

   mutationRate = 0.0;
   mutationRange = 0.0;

   // Stopping criteria

   evaluationGoal = -1.0e69;
   maximumTime = 1.0e69;

   maximumNumberOfGenerations = 100;

   showPeriod = 1;
}

// DESTRUCTOR

/// Destructor.

EvolutionaryAlgorithm::~EvolutionaryAlgorithm(void)
{

}


// METHODS

// int getPopulationSize(void) method

/// This method returns the number of individuals in the population.

int EvolutionaryAlgorithm::getPopulationSize(void)
{
   return(populationSize);
}


// void setMaximumNumberOfGenerations(int) method

/// This method sets a new maximum number of generations in the evolutionary 
/// algorithm optimization process.
///
/// @param newMaximumNumberOfGenerations:
/// Maximum number of generations.

void EvolutionaryAlgorithm
::setMaximumNumberOfGenerations(int newMaximumNumberOfGenerations)
{
   if(newMaximumNumberOfGenerations <= 0)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. "
                << "void setMaximumNumberOfGenerations(int) method."
                << std::endl
                << "Maximum number of generations must be greater than 0."
                << std::endl << std::endl;

      exit(1);
   }
   
   // Set maximum number of generations

   maximumNumberOfGenerations = newMaximumNumberOfGenerations;
}

// void setShowPeriod(int) method

/// This method sets a new number of generations between the optimization
/// showing progress.
///
/// @param newShowPeriod:
/// Number of generations between the optimization showing progress.

void EvolutionaryAlgorithm::setShowPeriod(int newShowPeriod)
{
   if(newShowPeriod <= 0)
   {
      std::cout << std::endl
                << "Error: void setShowPeriod(int) method."
                << std::endl
                << "Show period must be greater than 0."
                << std::endl << std::endl;

      exit(1);
   }
   
   // Set show period
   
   showPeriod = newShowPeriod;
}


// Vector<double> getMeanEvaluationHistory(void) method

/// This method returns a history with the mean evaluation of the population 
/// during optimization. 
///
/// @see getBestEvaluationHistory(void).
/// @see getStandardDeviationOfEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getMeanEvaluationHistory(void)
{
   return(meanEvaluationHistory);
}


// Vector<double> getStandardDeviationEvaluationHistory(void) method

/// This method returns a history with the standard deviation of the population  
/// evaluation during optimization. 
///
/// @see getBestEvaluationHistory(void).
/// @see getMeanEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getStandardDeviationEvaluationHistory(void)
{
   return(standardDeviationEvaluationHistory);
}


// Vector<double> getBestEvaluationHistory(void) method

/// This method returns a history with the evaluation value of the best individual 
/// ever during optimization. 
///
/// @see getMeanEvaluationHistory(void).
/// @see getStandardDeviationEvaluationHistory(void).

Vector<double> EvolutionaryAlgorithm::getBestEvaluationHistory(void)
{
   return(bestEvaluationHistory);
}


// int getShowPeriod(void) method

/// This method returns the number of generations between the 
/// evolutionary algorithm optimization showing progress. 

int EvolutionaryAlgorithm::getShowPeriod(void)
{
   return(showPeriod);    
}


// void setPopulationSize(int) method
//
/// This method sets a new population with a new number of individuals.  
/// The new population size must be an even number equal or greater than four. 
///
/// @param newPopulationSize: Number of individuals in the population.
/// This must be an even number equal or greater than four. 

void EvolutionaryAlgorithm::setPopulationSize(int newPopulationSize)
{
   if(newPopulationSize < 4)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. " << std::endl
                << "void setPopulationSize(int) method." << std::endl
                << "New population size must be equal or greater than 4."
                << std::endl << std::endl;

      exit(1);
   }
   else if(newPopulationSize%2 != 0)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. "
                << "void setPopulationSize(int) method." << std::endl
                << "New population size is not divisible by 2." << std::endl
                << std::endl;

      exit(1);
   }

   // Set population size

   populationSize = newPopulationSize;

   // Set population matrix

   Vector<double> lowerBound = objectiveFunction->getLowerBound();
   Vector<double> upperBound = objectiveFunction->getUpperBound();

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Matrix<double> newPopulation(populationSize, numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         newPopulation[i][j] = lowerBound[j] 
         + (upperBound[j] - lowerBound[j])*random;
      }
   }

   population = newPopulation;

   // Set evaluation vector

   Vector<double> newEvaluation(populationSize, 0.0);

   evaluation = newEvaluation;

   // Set fitness vector

   Vector<double> newFitness(populationSize, 0.0);

   fitness = newFitness;

   // Set selection vector

   Vector<bool> newSelection(populationSize, false);

   selection = newSelection;
}


// Population stuff

// Matrix<double> getPopulation(void) method

/// This method returns the population Matrix.

Matrix<double> EvolutionaryAlgorithm::getPopulation(void)
{
   return(population);
}


// Vector<double> getEvaluation(void) method

/// This method returns the actual evaluation value of all individuals in the 
/// population.
///
/// @see getFitness(void).
/// @see getSelection(void).

Vector<double> EvolutionaryAlgorithm::getEvaluation(void)
{
   return(evaluation);
}


// Vector<double> getFitness(void) method

/// This method returns the actual fitness value of all individuals in the population.
///
/// @see getEvaluation(void).
/// @see getSelection(void).

Vector<double> EvolutionaryAlgorithm::getFitness(void)
{
   return(fitness);
}


// Vector<bool> getSelection(void) method

/// This method returns the actual selection value of all individuals in the population.
///
/// @see getEvaluation(void).
/// @see getFitness(void).

Vector<bool> EvolutionaryAlgorithm::getSelection(void)
{
   return(selection);
}


// void setPopulation(Matrix<double>) method

/// This method sets a new population.
///
/// @param newPopulation: Population Matrix.

void EvolutionaryAlgorithm::setPopulation(Matrix<double> newPopulation)
{

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   if(newPopulation.getNumberOfRows() != populationSize)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. " << std::endl
                << "void setPopulation(Matrix<double>) method." << std::endl
                << "New population size is not equal to population size."
                << std::endl << std::endl;

      exit(1);
   }
   else if(newPopulation.getNumberOfColumns() != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. "
                << "void setPopulation(Matrix<double>) method." << std::endl
                << "New number of free parameters is not equal to "
                << "number of free parameters." << std::endl
                << std::endl;

      exit(1);
   }

   // Set population matrix

   population = newPopulation;
}


// void setEvaluation(Vector<double>) method

/// This method sets a new population evaluation Vector.
///
/// @param newEvaluation: Population evaluation values.

void EvolutionaryAlgorithm::setEvaluation(Vector<double> newEvaluation)
{
   if(newEvaluation.getSize() != populationSize)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. " << std::endl
                << "void setEvaluation(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl;

      exit(1);
   }

   // Set evaluation vector

   evaluation = newEvaluation;
}


// void setFitness(Vector<double>) method

/// This method sets a new population fitness Vector.
///
/// @param newFitness: Population fitness values.

void EvolutionaryAlgorithm::setFitness(Vector<double> newFitness)
{
   if(newFitness.getSize() != populationSize)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class. " << std::endl
                << "void setFitness(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl;

      exit(1);
   }

   // Set fitness vector

   fitness = newFitness;
}


// void setSelection(Vector<bool>) method

/// This method sets a new population selection Vector.
///
/// @param newSelection: Population selection values.

void EvolutionaryAlgorithm::setSelection(Vector<bool> newSelection)
{
   if(newSelection.getSize() != populationSize)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void setSelection(Vector<double>) method." << std::endl
                << "Size is not equal to population size." << std::endl
                << std::endl;

      exit(1);
   }
   
   // Set selection vector
   
   selection = newSelection;
}


// Vector<double> getIndividual(int) method

/// This method returns the Vector of free parameters corresponding to the 
/// individual i in the population.
///
/// @param i: Index of individual in the population.

Vector<double> EvolutionaryAlgorithm::getIndividual(int i)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   for(int j = 0; j < numberOfVariables; j++)
   {
      individual[j] = population[i][j];
   }

   return(individual);
}


// setIndividual(int, Vector<double>) method

/// This method sets a new Vector of free parameters to the individual i in the 
/// population. 
///
/// @param i: Index of individual in the population.
/// @param individual: Vector of free parameters to be assigned to individual i.

void EvolutionaryAlgorithm::setIndividual(int i, Vector<double> individual)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   for(int j = 0; j < numberOfVariables; j++)
   {
      population[i][j] = individual[j];
   }
}


// Optimization parameters

// double getSelectivePressure(void) method

/// This method returns the selective pressure value.
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see performNonLinearRankingFitnessAssignment(void).

double EvolutionaryAlgorithm::getSelectivePressure(void)
{
   return(selectivePressure);
}


// double getRecombinationSize(void) method

/// This method returns the recombination size value.
///
/// @see performIntermediateRecombination(void)
/// @see performLineRecombination(void).

double EvolutionaryAlgorithm::getRecombinationSize(void)
{
   return(recombinationSize);
}


// double getMutationRate(void) method

/// This method returns the mutation rate value.
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

double EvolutionaryAlgorithm::getMutationRate(void)
{
   return(mutationRate);
}


// double getMutationRange(void) method

/// This method returns the mutation range value.
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

double EvolutionaryAlgorithm::getMutationRange(void)
{
   return(mutationRange);
}


// FitnessAssignmentMethod getFitnessAssignmentMethod(void) method

/// This method returns the fitness assignment method used for optimization.
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see performNonLinearRankingFitnessAssignment(void).
 
EvolutionaryAlgorithm::FitnessAssignmentMethod
EvolutionaryAlgorithm::getFitnessAssignmentMethod(void)
{
   return(fitnessAssignmentMethod);
}


// SelectionMethod getSelectionMethod(void) method

/// This method returns the selection method used for optimization.
///
/// @see performRouletteWheelSelection(void).
/// @see performStochasticUniversalSamplingSelection(void).

EvolutionaryAlgorithm::SelectionMethod 
EvolutionaryAlgorithm::getSelectionMethod(void)
{
   return(selectionMethod);
}


// RecombinationMethod getRecombinationMethod(void) method

/// This method returns the recombination method used for optimization.
///
/// @see performIntermediateRecombination(void).
/// @see performLineRecombination(void).

EvolutionaryAlgorithm::RecombinationMethod EvolutionaryAlgorithm::getRecombinationMethod(void)
{
   return(recombinationMethod);
}


// MutationMethod getMutationMethod(void) method

/// This method returns the mutation method used for optimization.
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).

EvolutionaryAlgorithm::MutationMethod EvolutionaryAlgorithm::getMutationMethod(void)
{
   return(mutationMethod);
}


// void setSelectivePressure(double) method

/// This method sets a new value for the selective pressure parameter.
/// Linear ranking allows values for the selective pressure between 1 and 2.
/// Non linear ranking allows values for the selective pressure between 1
/// and [population size] - 2.
///
/// @param newSelectivePressure: Selective pressure value. This must be between 1 and 2 
/// for linear ranking fitness assignment and between 1 1 and [population size] - 2 
/// for non linear ranking fitness assignment. 
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see performNonLinearRankingFitnessAssignment(void).

void EvolutionaryAlgorithm::setSelectivePressure(double newSelectivePressure)
{
   switch(fitnessAssignmentMethod)
   {
      case LinearRanking:

         if((newSelectivePressure < 1.0) || (newSelectivePressure > 2.0))
         {
            std::cout << std::endl
                      << "Error: EvolutionaryAlgorithm class. " << std::endl
                      << "void setSelectivePressure(double) method. "
                      << "Case linear ranking." << std::endl
                      << "Selective pressure must be a value between 1 and 2."
                      << std::endl << std::endl;

            exit(1);
         }

         selectivePressure = newSelectivePressure;

      break;

      case NonLinearRanking:

         if((newSelectivePressure < 1.0) || (newSelectivePressure > populationSize-2))
         {
            std::cout << std::endl
                      << "Error: EvolutionaryAlgorithm class. "
                      << "void setSelectivePressure(double) method. "
                      << "Case non linear ranking." << std::endl
                      << "Selective pressure must be a value between 1 "
                      << "and population size - 2."
                      << std::endl << std::endl;

            exit(1);
         }

         selectivePressure = newSelectivePressure;

      break;
   }
}


// void setRecombinationSize(double) method

/// This method sets a new value for the recombination size parameter.
/// The recombination size value must be equal or greater than 0.
///
/// @param newRecombinationSize: Recombination size value. This must be equal or greater than 0.
///
/// @see performIntermediateRecombination(void)
/// @see performLineRecombination(void).

void EvolutionaryAlgorithm::setRecombinationSize(double newRecombinationSize)
{
   if(newRecombinationSize < 0.0)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void setRecombinationSize(double) method." << std::endl
                << "Recombination size must be equal or greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set recombination  size

   recombinationSize = newRecombinationSize;
}


// void setMutationRate(double) method

/// This method sets a new value for the mutation rate parameter.
/// The mutation rate value must be between 0 and 1.
///
/// @param newMutationRate: Mutation rate value. This must be equal or greater than 0
/// less or equal to 1. 
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

void EvolutionaryAlgorithm::setMutationRate(double newMutationRate)
{
   if(newMutationRate < 0.0 || newMutationRate > 1.0)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void setMutationRate(double) method." << std::endl
                << "Mutation rate must be a value between 0 and 1. " << std::endl
                << std::endl;

      exit(1);
   }

   // Set mutation rate

   mutationRate = newMutationRate;
}


// void setMutationRange(double) method

/// This method sets a new value for the mutation range parameter.
/// The mutation range value must be 0 or a positive number. 
///
/// @param newMutationRange: Mutation range value. This must be equal or greater than 0.
///
///  @see performNormalMutation(void).
///  @see performUniformMutation(void).

void EvolutionaryAlgorithm::setMutationRange(double newMutationRange)
{
   if(newMutationRange < 0.0)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void setMutationRange(double) method." << std::endl
                << "Mutation range must be equal or greater than 0. " << std::endl
                << std::endl;

      exit(1);
   }

   // Set mutation range

   mutationRange = newMutationRange;
}


// void setFitnessAssignmentMethod(FitnessAssignmentMethod) method

/// This method sets a new fitness assignment method to be used for optimization.
///
/// @param newFitnessAssignmentMethod: Fitness assignment method chosen for 
/// optimization.
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see performNonLinearRankingFitnessAssignment(void).

void EvolutionaryAlgorithm::setFitnessAssignmentMethod
(EvolutionaryAlgorithm::FitnessAssignmentMethod newFitnessAssignmentMethod)
{
   fitnessAssignmentMethod = newFitnessAssignmentMethod;
}


// void setSelectionMethod(SelectionMethod) method

/// This method sets a new selection method to be used for optimization.
///
/// @param newSelectionMethod: Selection method.
///
/// @see performRouletteWheelSelection(void).
/// @see performStochasticUniversalSamplingSelection(void).

void EvolutionaryAlgorithm::setSelectionMethod
(EvolutionaryAlgorithm::SelectionMethod newSelectionMethod)
{
   selectionMethod = newSelectionMethod;
}


// void setRecombinationMethod(RecombinationMethod) method

/// This method sets a new recombination method to be used for optimization.
///
/// @param newRecombinationMethod: Recombination method chosen for optimization. 
///
/// @see performIntermediateRecombination(void).
/// @see performLineRecombination(void).

void EvolutionaryAlgorithm::setRecombinationMethod
(EvolutionaryAlgorithm::RecombinationMethod newRecombinationMethod)
{
   recombinationMethod = newRecombinationMethod;
}


// void setMutationMethod(MutationMethod) method

/// This method sets a new mutation method to be used for optimization.
///
/// @param newMutationMethod: Mutation method chosen for optimization. 
///
/// @see performNormalMutation(void).
/// @see performUniformMutation(void).


void EvolutionaryAlgorithm::setMutationMethod
(EvolutionaryAlgorithm::MutationMethod newMutationMethod)
{
   mutationMethod = newMutationMethod;
}


// void initPopulationAtRandom(void) method
//
/// This method initializes all the individuals in the population at 
/// random, with values comprised between -1 and 1.
///
/// @see initPopulationAtRandom(double, double).
/// @see initPopulationAtRandom(Vector<double>, Vector<double>).
/// @see initPopulationAtRandom(Matrix<double>).

void EvolutionaryAlgorithm::initPopulationAtRandom(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = -1.0 + 2.0*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationAtRandom(double, double) method
//
/// This method initializes all the individuals in the population at 
/// random, with values comprised between a minimum and a maximum value.
///
/// @param minimumValue: Minimum initialization value.
/// @param maximumValue: Maximum initialization value.
///
/// @see initPopulationAtRandom(void).
/// @see initPopulationAtRandom(Vector<double>, Vector<double>).
/// @see initPopulationAtRandom(Matrix<double>).

void EvolutionaryAlgorithm
::initPopulationAtRandom(double minimumValue, double maximumValue)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = minimumValue + (maximumValue-minimumValue)*random;
      }

      setIndividual(i, individual);
   }
}


// void initPopulationAtRandom(Vector<double>, Vector<double>) method

/// This method initializes the free parameters of all the individuals in the population at 
/// random, with values comprised between different minimum and maximum values 
/// for each variable.
///
/// @param minimumValue: Minimum initialization value.
/// @param maximumValue: Maximum initialization value.
///
/// @see initPopulationAtRandom(void).
/// @see initPopulationAtRandom(double, double).
/// @see initPopulationAtRandom(Matrix<double>).

void EvolutionaryAlgorithm
::initPopulationAtRandom(Vector<double> minimumValue, Vector<double> maximumValue)
{
   int minimumValueSize = minimumValue.getSize();
   int maximumValueSize = maximumValue.getSize();

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   if(minimumValueSize != numberOfVariables || maximumValueSize != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void initPopulationAtRandom(Vector<double>, Vector<double>)."
                << std::endl
                << "Minimum value and maximum value sizes must be equal "
                << "to number of free parameters." 
                << std::endl << std::endl;
 
      exit(1);
   }

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = minimumValue[j] + (maximumValue[j]-minimumValue[j])*random;
      }

      setIndividual(i, individual);
   }
}

// void initPopulationAtRandom(Matrix<double>) method

/// This method initializes the free parameters of all the individuals in the population at 
/// random, with values comprised between different minimum and maximum values 
/// for each variable and from a single Matrix.
///
/// @param minimumAndMaximumValues: Minimum and maximum initialization values.
///
/// @see initPopulationAtRandom(void).
/// @see initPopulationAtRandom(double, double).
/// @see initPopulationAtRandom(Vector<double>, Vector<double>).

void EvolutionaryAlgorithm
::initPopulationAtRandom(Matrix<double> minimumAndMaximumValues)
{
   int numberOfRows = minimumAndMaximumValues.getNumberOfRows();
   int numberOfColumns = minimumAndMaximumValues.getNumberOfColumns();

   int numberOfVariables = objectiveFunction->getNumberOfVariables();


   if(numberOfRows != 2 || numberOfColumns != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class."  << std::endl
                << "void initPopulationAtRandom(Matrix<double>)." << std::endl
                << "Number of rows must be two and number of columns must be equal to "
                << "number of free parameters." 
                << std::endl << std::endl;
 
      exit(1);
   }

   Vector<double> minimumValue(numberOfVariables, 0.0);
   Vector<double> maximumValue(numberOfVariables, 0.0);

   for(int i = 0; i < numberOfVariables; i++)
   {
      minimumValue[i] = minimumAndMaximumValues[0][i];
      maximumValue[i] = minimumAndMaximumValues[1][i];
   }

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         double random = (double)rand()/(RAND_MAX+1.0);

         individual[j] = minimumValue[j] + (maximumValue[j]-minimumValue[j])*random;
      }

      setIndividual(i, individual);
   }
}


// void evaluatePopulation(void) method
//
/// This method evaluates the objective function of all individuals in the
/// population. Results are stored in the evaluation Vector.

void EvolutionaryAlgorithm::evaluatePopulation(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   // Evaluate the objective function for all individuals

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      evaluation[i] = objectiveFunction->getEvaluation(individual);
      
      if(!(evaluation[i] > -1.0e69 && evaluation[i] < 1.0e69))
      {
         std::cout << std::endl
                   << "Error: EvolutionaryAlgorithm class. " << std::endl
                   << "void evaluatePopulation(void) method." << std::endl
                   << "Evaluation of individual " << i 
                   << " is not a real number." << std::endl
                   << std::endl;

         exit(1);
      }          
   }
}


// void performLinearRankingFitnessAssignment(void) method
//
/// This method ranks all individuals in the population by their
/// evaluation value, so that the least fit individual has rank 1 and the
/// fittest individual has rank [population size].
/// It then assigns them a fitness value linearly proportional to their
/// rank. Results are stored in the fitness Vector.
///
/// @see performNonLinearRankingFitnessAssignment(void).

void EvolutionaryAlgorithm::performLinearRankingFitnessAssignment(void)
{
   // Sorted evaluation vector

   Vector<double> sortedEvaluation(populationSize, 0.0);

   sortedEvaluation = evaluation;

   std::sort(sortedEvaluation.begin(), sortedEvaluation.end(), std::less<double>());

   // Rank vector

   Vector<int> rank(populationSize, 0);

   for(int i = 0; i < populationSize; i++)
   {
      for(int j = 0; j < populationSize; j++)
      {
         if(evaluation[j] == sortedEvaluation[i])
         {
            rank[j] = populationSize - i;
         }
      }
   }

   // Perform linear ranking fitness assignment

   for(int i = 0; i < populationSize; i++)
   {
      fitness[i] = 2.0 - selectivePressure
      + 2.0*(selectivePressure - 1.0)*(rank[i] - 1.0)/(populationSize - 1.0);
      
      if(!(fitness[i] > -1.0e69 && fitness[i] < 1.0e69))
      {
         std::cout << std::endl
                   << "Error: EvolutionaryAlgorithm class. " << std::endl
                   << "void performLinearRankingFitnessAssignment(void) method." 
                   << std::endl
                   << "Fitness of individual " << i 
                   << " is not a real number." << std::endl
                   << std::endl;

         exit(1);
      }          
   }
}


// void performNonLinearRankingFitnessAssignment(void) method
//
/// This method ranks all individuals in the population by their 
/// evaluation value, so that the least fit individual has rank 1 and the 
/// fittest individual has rank population size.
/// Please, do not use this method. It is not yet implemented.
///
/// @see performLinearRankingFitnessAssignment(void).
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performNonLinearRankingFitnessAssignment(void)
{
   std::cout << std::endl
             << "Error: EvolutionaryAlgorithm class. " << std::endl
             << "void performNonLinearRankingFitnessAssignment(void) method." 
             << std::endl
             << "Sorry, this method is not yet implemented." << std::endl
             << std::endl;

   exit(1);
}


// void performRouletteWheelSelection(void) method

/// This metod performs selection with roulette wheel selection. It 
/// selects half of the individuals from the population. Results are
/// stored in the selection Vector. 
///
/// @see performStochasticUniversalSamplingSelection(void) method
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performRouletteWheelSelection(void)
{
   // Set selection vector to false 

   for(int i = 0; i < populationSize; i++)
   {
      selection[i] = false;
   }

   int numberOfSelectedIndividuals = populationSize/2;

   Vector<double> cumulativeFitness(populationSize, 0.0);

   // Cumulative fitness vector

   cumulativeFitness[0] = fitness[0]; 

   for(int i = 1; i < populationSize; i++)
   {
      cumulativeFitness[i] = cumulativeFitness[i-1] + fitness[i];
   }

   // Select individuals until the desired number of selections is obtained

   int countNumberOfSelectedIndividuals = 0;

   do
   {
      // Random number between 0 and total cumulative fitness

      double random = (double)rand()/(RAND_MAX+1.0);

      double pointer = cumulativeFitness[populationSize-1]*random;

      // Perform selection

      if(pointer < cumulativeFitness[0])
      {
         if(selection[0] == false)
         {
            selection[0] = true;
            countNumberOfSelectedIndividuals++;
         }
      }
      
      for(int i = 1; i < populationSize; i++)
      {
         if(pointer < cumulativeFitness[i] && pointer >= cumulativeFitness[i-1])
         {
            if(selection[i] == false)
            {
               selection[i] = true;
               countNumberOfSelectedIndividuals++;
            }
         }
      }
   }while(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals);

   // Control sentence

   if(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals)
   {
      std::cout << std::endl 
                << "Error: EvolutionaryAlgorithm class. " << std::endl
                << "void performRouletteWheelSelection(void) method."
                << std::endl
                << "Count number of selected individuals is not equal to "
                << "number of selected individuals."
                << std::endl << std::endl;

      exit(1);
   }
}


// void performStochasticUniversalSamplingSelection(void) method
//
/// This metod performs selection with stochastic universal sampling. It 
/// selects half of the individuals from the population. Results are
/// stored in the selection vector. 
///
/// @see performRouletteWheelSelection(void) method.
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performStochasticUniversalSamplingSelection(void)
{
   // Set selection vector to false

   for(int i = 0; i < populationSize; i++)
   {
      selection[i] = false;
   }
 
   int numberOfSelectedIndividuals = populationSize/2;

   Vector<double> cumulativeFitness(populationSize, 0.0);

   Vector<double> pointer(numberOfSelectedIndividuals, 0.0);

   // Cumulative fitness vector

   cumulativeFitness[0] = fitness[0];

   for(int i = 1; i < populationSize; i++)
   {  
      cumulativeFitness[i] = cumulativeFitness[i-1] + fitness[i];
   }


   // Pointer vector

   // Random number between 0 and 
   // totalCumulativeFitnees/numberOfSelectedIndividuals 

   double random = (double)rand()/(RAND_MAX+1.0);

   pointer[0] = random
   *cumulativeFitness[populationSize-1]/(double)numberOfSelectedIndividuals;

   for(int i = 1; i < numberOfSelectedIndividuals; i++)
   {
      pointer[i] = pointer[i-1] 
      + cumulativeFitness[populationSize-1]/(double)numberOfSelectedIndividuals;
   }

   // Selection vector

   int countNumberOfSelectedIndividuals = 0;

   if(pointer[0] <= cumulativeFitness[0])
   {
      selection[0] = true;
      countNumberOfSelectedIndividuals++;
   }


   for(int i = 0; i < numberOfSelectedIndividuals; i++)
   {

      for(int j = 1; j < populationSize; j++)
      {
         if(pointer[i] <= cumulativeFitness[j] && pointer[i] > cumulativeFitness[j-1])
         {
            selection[j] = true;
            countNumberOfSelectedIndividuals++;
         }
      }
   }

   // Number of selected individuals control sentence 

   if(countNumberOfSelectedIndividuals != numberOfSelectedIndividuals)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void performStochasticUniversalSamplingSelection(void) method."
                << std::endl
                << "Count number of selected individuals is not equal to "
                << "number of selected individuals." 
                << std::endl << std::endl;

      exit(1);
   }
}


// void performIntermediateRecombination(void) method

/// This method performs intermediate recombination between pairs
/// of selected individuals to generate a new population. Each selected
/// individual is to be recombined with two other selected individuals
/// chosen at random. Results are stored in the population matrix.
///
/// @see recombinationSize.
///
/// @see performLineRecombination(void).
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performIntermediateRecombination(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Matrix<double> newPopulation(populationSize, numberOfVariables, 0.0);

   int numberOfSelectedIndividuals = populationSize/2;

   int count = 0;
   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
         count ++;        
   }

   Vector<double> parent1(numberOfVariables, 0.0);
   Vector<double> parent2(numberOfVariables, 0.0);

   Vector<double> offspring(numberOfVariables, 0.0);

   Matrix<int> recombination(populationSize, 2, 0);

   // Start recombination	

   int countNewPopulationSize = 0;

   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
      {
         // Set parent 1

         parent1 = getIndividual(i);

         // Generate 2 offspring with parent 1

         for(int j = 0; j < 2; j++)
         {
            // Choose parent 2 at random among selected individuals	

            bool parent2Candidate = false;

            do{
               // Integer random number beteen 0 and population size

               double random = (double)rand()/(RAND_MAX+1.0);

               int parent2CandidateIndex = (int)(populationSize*random);

               // Check if candidate for parent 2 is ok

               if(selection[parent2CandidateIndex] == true 
               && parent2CandidateIndex != i)
               {
                  parent2Candidate = true;

                  recombination[countNewPopulationSize][0] = i;

                  recombination[countNewPopulationSize][1]
                  = parent2CandidateIndex;

                  parent2 = getIndividual(parent2CandidateIndex);

                  // Perform intermediate recombination between parent 1
                  // and parent 2

                  for(int j = 0; j < numberOfVariables; j++)
                  {
                     // Choose the scaling factor to be a random number between
                     // -recombinationSize and 1+recombinationSize for each
                     // variable anew.

                     double random = (double)rand()/(RAND_MAX+1.0);

                     double scalingFactor = -1.0*recombinationSize
                     + (1.0 + recombinationSize)*random;

                     offspring[j] = scalingFactor*parent1[j]
                     + (1.0 - scalingFactor)*parent2[j];
                  }

                  // Add offspring to newPopulation matrix

                  for(int j = 0; j < numberOfVariables; j++)
                  {
                     newPopulation[countNewPopulationSize][j] = offspring[j];	
                  }

                  countNewPopulationSize++;
               }
            }while(parent2Candidate != true);
         }
      }
   }

   // Count number of new individuals control sentence

   if(countNewPopulationSize != populationSize)
   {
      std::cout << std::endl 
                << "Error: EvolutionaryAlgorithm class. "
                << "void performIntermediateRecombination(void) method." 
                << std::endl
                << "Count new population size is not equal to population size." 
                << std::endl
                << std::endl;

      exit(1);
   }

   // Set new population

   population = newPopulation;
}


// void performLineRecombination(void) method

/// This method performs line recombination between pairs
/// of selected individuals to generate a new population. Each selected
/// individual is to be recombined with two other selected individuals
/// chosen at random. Results are stored in the population matrix.
///
/// @see recombinationSize.
///
/// @see performIntermediateRecombination(void).
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performLineRecombination(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Matrix<double> newPopulation(populationSize, numberOfVariables, 0.0);

   int numberOfSelectedIndividuals = populationSize/2;

   Vector<double> parent1(numberOfVariables, 0.0);
   Vector<double> parent2(numberOfVariables, 0.0);

   Vector<double> offspring(numberOfVariables, 0.0);

   Matrix<int> recombination(populationSize, 2, 0);

   // Start recombination	

   int countNewPopulationSize = 0;

   for(int i = 0; i < populationSize; i++)
   {
      if(selection[i] == true)
      {

         // Set parent 1

         parent1 = getIndividual(i);

         // Generate 2 offspring with parent 1

         for(int j = 0; j < 2; j++)
         {
            // Choose parent 2 at random among selected individuals	

            bool parent2Candidate = false;

            do
            {
               // Integer random number beteen 0 and population size

               double random = (double)rand()/(RAND_MAX + 1.0);

               int parent2CandidateIndex = (int)(populationSize*random);

               // Check if candidate for parent 2 is ok

               if(selection[parent2CandidateIndex] == true && parent2CandidateIndex != i)
               {
                  parent2Candidate = true;

                  recombination[countNewPopulationSize][0] = i;
                  recombination[countNewPopulationSize][1] = parent2CandidateIndex;

                  parent2 = getIndividual(parent2CandidateIndex);

                  // Perform intermediate recombination between parent 1
                  // and parent 2

                  // Choose the scaling factor to be a random number between
                  // -recombinationSize and 1+recombinationSize for all
                  // variables.

                  double random = (double)rand()/(RAND_MAX+1.0);

                  double scalingFactor = -1.0*recombinationSize 
                  + (1.0 + recombinationSize)*random;

                  for(int j = 0; j < numberOfVariables; j++)
                  {
                     offspring[j] = scalingFactor*parent1[j]
                     + (1.0 - scalingFactor)*parent2[j];
                  }

                  // Add offspring to newPopulation matrix

                  for(int j = 0; j < numberOfVariables; j++)
                  {
                     newPopulation[countNewPopulationSize][j] = offspring[j];	
                  }

                  countNewPopulationSize++;
               }
            }while(parent2Candidate == false);
         }
      }
   }

   // Count new population size control sentence

   if(countNewPopulationSize != populationSize)
   {
      std::cout << std::endl
                << "Error:  EvolutionaryAlgorithm class. " << std::endl
                << "void performLineRecombination(void) method."
                << std::endl
                << "Count new population size is not equal to population size."
                << std::endl << std::endl;

      exit(1);
   }

   // Set new population

   population = newPopulation;
}


// void performNormalMutation(void) method

/// This method performs normal mutation to all individuals in order to generate
/// a new population. It uses the Box-Muller transformation to generate
/// random numbers with normal distribution. Results are stored in the
/// population matrix.
///
/// @see mutationRate.
/// @see mutationRange.
///
/// @see performUniformMutation(void).
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performNormalMutation(void)
{
   const double pi = 3.141592654;

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         // Random number between 0 and 1

         double pointer = (double)rand()/(RAND_MAX+1.0);

         if(pointer < mutationRate)
         {
            // Random numbers between 0 and 1

            double random1 = 0.0;
            
            do // random1 cannot be zero
            {
               random1 = (double)rand()/(RAND_MAX+1.0);
            
            }while(random1 == 0);

            double random2 = (double)rand()/(RAND_MAX+1.0);

            // Box-Muller transformation

            double mean = 0.0;
            double standardDeviation = mutationRange;

            double normallyDistributedRandomNumber 
            = mean + sqrt(-2.0*log(random1))*sin(2.0*pi*random2)*standardDeviation;

            individual[j] += normallyDistributedRandomNumber;
         }
      }

      setIndividual(i, individual);
   }
}  


// void performUniformMutation(void) method

/// This method performs uniform mutation to all individuals in order to generate
/// a new population. Results are stored in the population matrix.
///
/// @see mutationRate.
/// @see mutationRange.
///
/// @see performNormalMutation(void).
/// @see getMinimalArgument(void).

void EvolutionaryAlgorithm::performUniformMutation(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> individual(numberOfVariables, 0.0);

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      for(int j = 0; j < numberOfVariables; j++)
      {
         // random number between 0 and 1

         double pointer = (double)rand()/(RAND_MAX+1.0);

         if(pointer < mutationRate)
         {
            // random number between 0 and 1

            double random = (double)rand()/(RAND_MAX+1.0);

            double uniformlyDistributedRandomNumber
            = (-1.0 + 2.0*random)*mutationRange;

            individual[j] += uniformlyDistributedRandomNumber;
         }
      }

      setIndividual(i, individual);
   }
}


// Vector<double> EvolutionaryAlgorithm::getMinimalArgument(void) method

/// This method optimizes the objective function with the evolutionary algorithm.
/// This process occurs according to the optimization operators and their related
/// parameters.
///
/// @see FitnessAssignmentMethod.
/// @see SelectionMethod.
/// @see RecombinationMethod.
/// @see MutationMethod.

Vector<double> EvolutionaryAlgorithm::getMinimalArgument(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> minimalArgument(numberOfVariables, 0.0);

   // Mean evaluation history vector

   Vector<double> newMeanEvaluationHistory(maximumNumberOfGenerations+1, 0.0);

   meanEvaluationHistory = newMeanEvaluationHistory;

   // Standard deviation of evaluation history vector

   Vector<double> newStandardDeviationEvaluationHistory(maximumNumberOfGenerations+1, 0.0);

   standardDeviationEvaluationHistory = newStandardDeviationEvaluationHistory;

   // Best evaluation history vector

   Vector<double> newBestEvaluationHistory(maximumNumberOfGenerations+1, 0.0);

   bestEvaluationHistory = newBestEvaluationHistory;

   Vector<double> bestIndividual(numberOfVariables, 0.0);

   double bestEvaluation = 1.0e69;

   time_t beginningTime, currentTime;

   double elapsedTime = 0.0;

   // Set beginning time

   time(&beginningTime);

   std::cout << std::endl
             << "Getting minimal argument with evolutionary algorithm..." 
             << std::endl;

   // Evaluation of population

   evaluatePopulation();

   // Check for best individual

   for(int i = 0; i < populationSize; i++)
   {
      if(evaluation[i] < bestEvaluation)
      {
         bestEvaluation = evaluation[i];

         bestIndividual = getIndividual(i);
      }
   }

   bestEvaluationHistory[0] = bestEvaluation;

   // Calculate mean population evaluation

   double meanEvaluation = 0.0;

   double sum1 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      sum1 += evaluation[i];
   }

   meanEvaluation = sum1/(double)populationSize;

   meanEvaluationHistory[0] = meanEvaluation;

   // Calculate standard deviation of evaluation

   double standardDeviationEvaluation = 0.0;

   double sum2 = 0.0;

   for(int i = 0; i < populationSize; i++)
   {
      sum2 += pow(evaluation[i] - meanEvaluation, 2);
   }

   standardDeviationEvaluation = sqrt(sum2/(double)populationSize);

   standardDeviationEvaluationHistory[0] = standardDeviationEvaluation;

   if(bestEvaluation <= evaluationGoal)
   {          
      std::cout << std::endl
                << "Initial evaluation is less than goal." << std::endl
                << "Initial evaluation: " << bestEvaluation << std::endl;

      objectiveFunction->print(bestIndividual);
      
      minimalArgument = bestIndividual;

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
      std::cout << std::endl
                   << "Initial Population:" << std::endl
                   << "Mean evaluation: " << meanEvaluation << std::endl
                   << "Standard deviation of evaluation: " << standardDeviationEvaluation
                   << std::endl
                   << "Best individual: " << std::endl
                   << "Evaluation: " << bestEvaluation << ";" << std::endl;

      objectiveFunction->print(bestIndividual);
   }

   // Main loop

   for(int generation = 1; generation <= maximumNumberOfGenerations; generation++)
   {
      // Fitness assignment

      switch(fitnessAssignmentMethod)
      {
         case LinearRanking:

            performLinearRankingFitnessAssignment();

            break;

         case NonLinearRanking:

            performNonLinearRankingFitnessAssignment();

            break;
      }

      // Selection

      switch(selectionMethod)
      {
         case RouletteWheel:

            performRouletteWheelSelection();

            break;

         case StochasticUniversalSampling:

            performStochasticUniversalSamplingSelection();

            break;
      }

      // Recombination

      switch(recombinationMethod)
      {
         case Intermediate:

            performIntermediateRecombination();

            break;

         case Line:

            performLineRecombination();

            break;
      }

      // Mutation

      switch(mutationMethod)
      {
         case Normal:

            performNormalMutation();

            break;

         case Uniform:

            performUniformMutation();

            break;
      }


      // Evaluation of population

      evaluatePopulation();

      // Check for best individual

      for(int i = 0; i < populationSize; i++)
      {
         if(evaluation[i] < bestEvaluation)
         {
            bestEvaluation = evaluation[i];

            bestIndividual = getIndividual(i);
         }
      }

      bestEvaluationHistory[generation] = bestEvaluation;

      // Calculate mean population evaluation

      double meanEvaluation = 0.0;

      double sum1 = 0.0;

      for(int i = 0; i < populationSize; i++)
      {
         sum1 += evaluation[i];
      }

      meanEvaluation = sum1/(double)populationSize;

      meanEvaluationHistory[generation] = meanEvaluation;

      // Calculate standard deviation of evaluation

      double standardDeviationEvaluation = 0.0;

      double sum2 = 0.0;

      for(int i = 0; i < populationSize; i++)
      {
         sum2 += pow(evaluation[i] - meanEvaluation, 2);
      }

      standardDeviationEvaluation = sqrt(sum2/(double)populationSize);

      standardDeviationEvaluationHistory[generation] = standardDeviationEvaluation;

      // Stopping Criteria

      // Evaluation goal 

      if (bestEvaluation <= evaluationGoal)
      {
         std::cout << std::endl
                   << "Generation " << generation << ": "
                   << "Evaluation goal reached." << std::endl;

         std::cout << "Final evaluation:" << std::endl;

         std::cout << "Evaluation: " << bestEvaluation << ";" << std::endl;
         
         objectiveFunction->print(bestIndividual);

         break;
      }

      // Maximum time

      time(&currentTime);

      elapsedTime = difftime(currentTime, beginningTime);

      if (elapsedTime >= maximumTime)
      {
         std::cout << std::endl
                   << "Generation " << generation << ": "
                   << "Maximum time reached."
                   << std::endl;

         std::cout << "Final evaluation:" << std::endl;

         std::cout << "Evaluation: " << bestEvaluation << ";" << std::endl;

         objectiveFunction->print(bestIndividual);

         break;
      }

      // Maximum number of generations

      if (generation == maximumNumberOfGenerations)
      {
         std::cout << std::endl
                   << "Generation " << generation << ": "
                   << "Maximum number of generations reached."
                   << std::endl;

         std::cout << "Final evaluation:" << std::endl;

         std::cout << "Evaluation: " << bestEvaluation << ";" << std::endl;

         objectiveFunction->print(bestIndividual);

         break;
      }

      // Progress

      if(generation % showPeriod == 0)
      {
         std::cout << std::endl
                   << "Generation " << generation << "; " << std::endl
                   << "Population: " << std::endl
                   << "Mean evaluation: " << meanEvaluation << std::endl
                   << "Standard deviation of evaluation: " << standardDeviationEvaluation
                   << std::endl
                   << "Best individual: " << std::endl;

         std::cout << "Evaluation: " << bestEvaluation << ";" << std::endl;
         
         objectiveFunction->print(bestIndividual);
      }


      // Reset selection vector

      Vector<bool> newSelection(populationSize, false);

      selection = newSelection;
   }

   // Set minimal argument

   minimalArgument = bestIndividual;

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

/// This method prints to the screen the members of the evolutionary algorithm 
/// object.
///
/// Optimization operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations.
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Generations between showing progress.
/// </ul>

void EvolutionaryAlgorithm::print(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   std::cout << std::endl
             << "Flood Neural Network. Evolutionary Algorithm Object." 
             << std::endl;

   std::cout << std::endl
             << "Population size: " << std::endl
             << populationSize << std::endl
             << "Number of free parameters: " << std::endl
             << numberOfVariables  << std::endl;

   // Optimization operators

   std::cout << std::endl
             << "Optimization operators:" << std::endl;

   // Fitness assingment method

   std::cout << "Fitness assignment method:" << std::endl;

   if(fitnessAssignmentMethod == LinearRanking)
   {
      std::cout << "Linear ranking" << std::endl;
   }
   else if(fitnessAssignmentMethod == NonLinearRanking)
   {
      std::cout << "Non linear ranking" << std::endl;
   }

   // Selection method

   std::cout << "Selection method:" << std::endl;

   if(selectionMethod == RouletteWheel)
   {
      std::cout << "Roulette wheel" << std::endl;
   }
   else if(selectionMethod == StochasticUniversalSampling)
   {
      std::cout << "Stochastic universal sampling" << std::endl;
   }

   // Recombination method

   std::cout << "Recombination method:" << std::endl;

   if(recombinationMethod == Line)
   {
      std::cout << "Line" << std::endl;
   }
   else if(recombinationMethod == Intermediate)
   {
      std::cout << "Intermediate" << std::endl;
   }

   // Mutation method

   std::cout << "Mutation method:" << std::endl;

   if(mutationMethod == Normal)
   {
      std::cout << "Normal" << std::endl;
   }
   else if(mutationMethod == Uniform)
   {
      std::cout << "Uniform" << std::endl;
   }


   // Optimization parameters

   std::cout << std::endl
             << "Optimization parameters: " << std::endl
             << "Selective pressure: " << std::endl
             << selectivePressure << std::endl
             << "Recombination size: " << std::endl
             << recombinationSize << std::endl
             << "Mutation rate: " << std::endl
             << mutationRate << std::endl
             << "Mutation range: " << std::endl
             << mutationRange << std::endl;

   // Stopping criteria

   std::cout << std::endl
             << "Stopping criteria: " << std::endl
             << "Evaluation goal: " << std::endl
             << evaluationGoal << std::endl
             << "Maximum time: " << std::endl
             << maximumTime << std::endl
             << "Maximum number of generations: " << std::endl
             << maximumNumberOfGenerations << std::endl;

   // User stuff

   std::cout << std::endl
             << "User stuff: " << std::endl
             << "Show period: " << showPeriod
             << std::endl;

   Vector<double> individual(numberOfVariables, 0.0);

   std::cout << std::endl
             << "Population:" << std::endl;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      std::cout << "Individual " << i << ":" << std::endl; 

      for(int j = 0; j < numberOfVariables; j++)
      {
         std::cout << individual[j] << " ";
      }

      std::cout << std::endl;
   }
}


// void save(char*) method

/// This method saves the evolutionary algorithm object to a data file. 
///
/// Optimization operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations.
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Generations between showing progress.
/// </ul>
///
/// @param filename: Filename.
///
/// @see load(char*).

void EvolutionaryAlgorithm::save(char* filename)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void save(char*) method." 
                << std::endl
                << "Cannot open evolutionary algorithm object data file."  
                << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving evolutionary algorithm object to data file..." 
                << std::endl;
   }

   // Write file header

   // Write file header

   file << "% Purple: An Open Source Numerical Optimization C++ Library." 
        << std::endl 
        << "% Newton Method Object." << std::endl; 

   file << "PopulationSize:" << std::endl
        << populationSize << std::endl
        << "NumberOfVariables:" << std::endl
        << numberOfVariables  << std::endl;

   // Optimization operators

   // Fitness assingment method

   file << "FitnessAssignmentMethod:" << std::endl;

   if(fitnessAssignmentMethod == LinearRanking)
   {
      file << "LinearRanking" << std::endl;
   }
   else if(fitnessAssignmentMethod == NonLinearRanking)
   {
      file << "NonLinearRanking" << std::endl;
   }

   // Selection method

   file << "SelectionMethod:" << std::endl;

   if(selectionMethod == RouletteWheel)
   {
      file << "RouletteWheel" << std::endl;
   }
   else if(selectionMethod == StochasticUniversalSampling)
   {
      file << "StochasticUniversalSampling" << std::endl;
   }

   // Recombination method

   file << "RecombinationMethod:" << std::endl;

   if(recombinationMethod == Line)
   {
      file << "Line" << std::endl;
   }
   else if(recombinationMethod == Intermediate)
   {
      file << "Intermediate" << std::endl;
   }

   // Mutation method

   file << "MutationMethod:" << std::endl;

   if(mutationMethod == Normal)
   {
      file << "Normal" << std::endl;
   }
   else if(mutationMethod == Uniform)
   {
      file << "Uniform" << std::endl;
   }

   // Optimization parameters

   file << "SelectivePressure:" << std::endl
        << selectivePressure << std::endl
        << "RecombinationSize:" << std::endl
        << recombinationSize << std::endl
        << "MutationRate:" << std::endl
        << mutationRate << std::endl
        << "MutationRange: " << std::endl
        << mutationRange << std::endl;

   // Stopping criteria

   file << "EvaluationGoal: " << std::endl
        << evaluationGoal << std::endl
        << "MaximumTime: " << std::endl
        << maximumTime << std::endl
        << "MaximumNumberOfGenerations: " << std::endl
        << maximumNumberOfGenerations << std::endl;

   // User stuff

   file << "ShowPeriod: " << std::endl
        << showPeriod << std::endl;

   Vector<double> individual(numberOfVariables, 0.0);

   file << "Population:" << std::endl;

   for(int i = 0; i < populationSize; i++)
   {
      individual = getIndividual(i);

      file << "Individual" << i << ":" << std::endl; 

      for(int j = 0; j < numberOfVariables; j++)
      {
         file << individual[j] << " ";
      }

      file << std::endl;
   }

   file.close();
}


// void load(char*) method

/// This method loads an evolutionary $algorithm object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
/// Optimization operators:
/// <ul>
/// <li> Fitness assignment method.
/// <li> Selection method.
/// <li> Recombination method.
/// <li> Mutation method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> Population size.
/// <li> Selective pressure.
/// <li> Recombination size.
/// <li> Mutation rate.
/// <li> Mutation range.
/// </ul>
///
/// Stopping criteria:
/// <ul> 
/// <li> Evaluation goal.
/// <li> Maximum time.
/// <li> Maximum number of generations.
/// </ul> 
///  
/// User stuff: 
/// <ul>
/// <li> Generations between showing progress.
/// </ul>
///
/// @param filename: Filename.
///
/// @see save(char*).

void EvolutionaryAlgorithm::load(char* filename)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." 
                << std::endl
                << "Cannot open evolutionary algorithm object data file."  
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading evolutionary algorithm object from data file..."  
                << std::endl;
   }

   // Load new number of individuals and new number of free parameters form
   // file

   int newPopulationSize = 0;
   int newNumberOfVariables = 0;

   std::string word;

   // Optimization parameters

   // Population size

   while(word != "PopulationSize:")
   {
      file >> word;
   }

   file >> newPopulationSize;

   if(newPopulationSize != populationSize)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "New population size is not equal to population size." << std::endl
                << std::endl;

      exit(1);
   }

   // Number of free parameters

   file >> word;

   file >> newNumberOfVariables;

   if(newNumberOfVariables != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: EvolutionaryAlgorithm class." << std::endl
                << "void load(char*) method." << std::endl
                << "New number of free parameters is not equal to number of free parameters."
                << std::endl << std::endl;

      exit(1);
   }


   // Optimization operators

   // Fitness assingment method

   file >> word;

   file >> word;

   if(word == "LinearRanking")
   {
      fitnessAssignmentMethod = LinearRanking;
   }
   else if(word == "NonLinearRanking")
   {
      fitnessAssignmentMethod = NonLinearRanking;
   }
 
   // Selection method
 
   file >> word;

   file >> word;

   if(word == "RouletteWheel")
   {
      selectionMethod = RouletteWheel;
   }
   else if(word == "StochasticUniversalSampling")
   {
      selectionMethod = StochasticUniversalSampling;
   } 

   // Recombination method

   file >> word;

   file >> word;

   if(word == "Line")
   {
      recombinationMethod = Line;
   }
   else if(word == "Intermediate")
   {
      recombinationMethod = Intermediate;
   } 

   // Mutation method

   file >> word;

   file >> word;

   if(word == "Normal")
   {
      mutationMethod = Normal;
   }
   else if(word == "Uniform")
   {
      mutationMethod = Uniform;
   } 

   // Optimization parameters

   // Selective pressure

   file >> word;

   file >> selectivePressure;

   // Recombination size

   file >> word;

   file >> recombinationSize;

   // Mutation rate

   file >> word;

   file >> mutationRate;

   // Mutation range
 
   file >> word;

   file >> mutationRange;

   // Stopping criteria: 

   // Evaluation goal

   file >> word;

   file >> evaluationGoal;

   // Maximum time

   file >> word;

   file >> maximumTime;
   
   // Maximum number of generations

   file >> word;

   file >> maximumNumberOfGenerations;

   // User stuff: 

   // Show period

   file >> word;

   file >> showPeriod;

   // Population

   file >> word;

   for (int i = 0; i < populationSize; i++)
   {
      file >> word;

      for (int j = 0; j < numberOfVariables; j++)
      {
         file >> population[i][j];
      }
   }

   // Close file

   file.close();
}


// void saveOptimizationHistory(char*) method

/// This method saves the optimization history to a data file:
///
/// <ol>
/// <li> Generation.
/// <li> Mean evaluation.
/// <li> Standard deviation of evaluation.
/// <li> Best evaluation ever.
/// </ol>
///
/// @param filename: Filename.

void EvolutionaryAlgorithm::saveOptimizationHistory(char* filename)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   int numberOfGenerations = meanEvaluationHistory.getSize();

   std::fstream file; 

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: EvolutionaryAlgorithm class. " << std::endl
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
        << "% Evolutionary Algorithm Optimization History." << std::endl
        << "% 1 - Generation." << std::endl
        << "% 2 - Standard Deviation of Evaluation." << std::endl
        << "% 3 - Best Evaluation Ever." << std::endl;

   // Write file data

   int size = bestEvaluationHistory.getSize();

   for (int i = 0; i < size; i++)
   {
      file << i << " "
           << meanEvaluationHistory[i] << " "
           << standardDeviationEvaluationHistory[i] << " "
           << bestEvaluationHistory[i] << std::endl;
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
