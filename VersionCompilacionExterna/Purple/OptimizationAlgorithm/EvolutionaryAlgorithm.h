/******************************************************************************/
/*                                                                            */
/*   E V O L U T I O N A R Y   A L G O R I T H M   C L A S S   H E A D E R    */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#ifndef __EVOLUTIONARYALGORITHM_H__
#define __EVOLUTIONARYALGORITHM_H__


#include "OptimizationAlgorithm.h"
#include "../ObjectiveFunction/ObjectiveFunction.h"

namespace Purple
{
 
/// This concrete class represents an evolutionary optimization algorithm
/// for an objective function.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

class EvolutionaryAlgorithm : public OptimizationAlgorithm
{

public:

   // ENUMERATIONS

   /// Enumeration of the available optimization operators for fitness assignment.

   enum FitnessAssignmentMethod{LinearRanking, NonLinearRanking};

   /// Enumeration of the available optimization operators for selection. 

   enum SelectionMethod{RouletteWheel, StochasticUniversalSampling};

   /// Enumeration of the available optimization operators for recombination.

   enum RecombinationMethod{Line, Intermediate};

   /// Enumeration of the available optimization operators for mutation.

   enum MutationMethod{Normal, Uniform};

private:

   // FIELDS

   // Population stuff

   /// Number of individuals in the population. 

   int populationSize;

   /// Population matrix.

   Matrix<double> population;

   /// Function evaluation of population.

   Vector<double> evaluation;

   /// Fitness of population.

   Vector<double> fitness;

   /// Selected individuals in population.

   Vector<bool> selection;

   // Optimization parameters

   /// Selective pressure. 
   /// Linear ranking allows values for the selective pressure between 1 and 2.
   /// Non linear ranking allows values for the selective pressure 
   /// between 1 and [population size] - 2.

   double selectivePressure;

   /// Recombination size. 
   /// The recombination size value must be equal or greater than 0.

   double recombinationSize;

   /// Mutation rate.
   /// The mutation rate value must be between 0 and 1.

   double mutationRate;

   /// Mutation range.
   /// The mutation range value must be 0 or a positive number. 

   double mutationRange;

   // Stopping criteria

   /// Maximum number of generations.

   int maximumNumberOfGenerations;

   /// Number of generations between showing progress.

   int showPeriod;

   /// Mean evaluation optimization history.

   Vector<double> meanEvaluationHistory;

   /// Standard deviation of evaluation optimization history.

   Vector<double> standardDeviationEvaluationHistory;

   /// Best evaluation ever optimization history.

   Vector<double> bestEvaluationHistory;


   /// Fitness assignment optimization operators enumeration.

   FitnessAssignmentMethod fitnessAssignmentMethod;

   /// Selection optimization operators enumeration.

   SelectionMethod selectionMethod;

   /// Recombination optimization operators enumeration.

   RecombinationMethod recombinationMethod;

   /// Mutation optimization operators enumeration.

   MutationMethod mutationMethod;

public:

   // GENERAL CONSTRUCTOR

   EvolutionaryAlgorithm(ObjectiveFunction*);


   // DEFAULT CONSTRUCTOR

   EvolutionaryAlgorithm(void);


   // DESTRUCTOR

   virtual ~EvolutionaryAlgorithm(void);


   // METHODS

   // Get methods

   int getPopulationSize(void);

   Matrix<double> getPopulation(void);

   Vector<double> getEvaluation(void);
   Vector<double> getFitness(void);
   Vector<bool> getSelection(void);

   double getSelectivePressure(void);
   double getRecombinationSize(void);
   double getMutationRate(void);
   double getMutationRange(void);

   int getMaximumNumberOfGenerations(void);
   int getShowPeriod(void);

   Vector<double> getMeanEvaluationHistory(void);
   Vector<double> getStandardDeviationEvaluationHistory(void);
   Vector<double> getBestEvaluationHistory(void);

   FitnessAssignmentMethod getFitnessAssignmentMethod(void);
   SelectionMethod getSelectionMethod(void);
   RecombinationMethod getRecombinationMethod(void);
   MutationMethod getMutationMethod(void);

   // Set methods

   void setPopulationSize(int);

   void setPopulation(Matrix<double>);

   void setEvaluation(Vector<double>);
   void setFitness(Vector<double>);
   void setSelection(Vector<bool>);

   void setSelectivePressure(double);
   void setRecombinationSize(double);

   void setMutationRate(double);
   void setMutationRange(double);

   void setFitnessAssignmentMethod(FitnessAssignmentMethod);
   void setSelectionMethod(SelectionMethod);
   void setRecombinationMethod(RecombinationMethod);
   void setMutationMethod(MutationMethod);

   void setMaximumNumberOfGenerations(int);
   void setShowPeriod(int);

   // Population methods

   Vector<double> getIndividual(int);
   void setIndividual(int, Vector<double>);

   void initPopulationAtRandom(void);
   void initPopulationAtRandom(double, double);
   void initPopulationAtRandom(Vector<double>, Vector<double>);
   void initPopulationAtRandom(Matrix<double>);

   // Population evaluation methods

   void evaluatePopulation(void);

   // Fitness assignment methods

   void performLinearRankingFitnessAssignment(void);
   void performNonLinearRankingFitnessAssignment(void);

   // Selection methods

   void performRouletteWheelSelection(void);
   void performStochasticUniversalSamplingSelection(void);


   // Recombination methods

   void performIntermediateRecombination(void);
   void performLineRecombination(void);

   // Mutation methods

   void performNormalMutation(void);
   void performUniformMutation(void);

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
