/******************************************************************************/
/*                                                                            */
/*   C O N J U G A T E   G R A D I E N T   C L A S S                          */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/

#include "ConjugateGradient.h"

#include <iostream>
#include <fstream>
#include <algorithm>
#include <functional>
#include <limits>
#include <math.h>
#include <time.h>

namespace Purple
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates a conjugate gradient optimization
/// algorithm object associated to an objective function object. 
/// It also initializes the class members to their default values:
///
/// Initial argument: Random point whithin the objective function domain.
///
/// Optimization operators:
/// <ul>
/// <li> Search direction method = Polak-Ribiere;
/// <li> Optimal step size method = Brent method;
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> First step size: 1.0e-3.
/// <li> Optimal step size tolerance: 1.0e-3.
/// </ul>
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
/// <li> Warning step size: 1000.
/// <li> Show period: 25. 
/// </ul>
///
/// @param newObjectiveFunction: Pointer to an objective function object.
///
/// @see ObjectiveFunction.
/// @see OptimizationAlgorithm.

ConjugateGradient::ConjugateGradient(ObjectiveFunction* newObjectiveFunction)
: OptimizationAlgorithm(newObjectiveFunction)
{
   // Optimization operators

   searchDirectionMethod = PolakRibiere;
   optimalStepSizeMethod = BrentMethod;

   // Initial argument

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

   // Optimization parameters 

   firstStepSize = 1.0e-3;

   optimalStepSizeTolerance = 1.0e-3;

   // Stopping criteria

   evaluationGoal = -1.0e69;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 1000;

   // User stuff

   warningStepSize = 1000.0;
   showPeriod = 25;
}


// DEFAULT CONSTRUCTOR

/// Default constructor. It creates a conjugate gradient optimization algorithm 
/// object not associated to any objective function object. 
/// It also initializes the class members to their default values:
///
/// Optimization operators:
/// <ul>
/// <li> Search direction method = Polak-Ribiere;
/// <li> Optimal step size method = Brent method;
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> First step size: 1.0e-3.
/// <li> Optimal step size tolerance: 1.0e-3.
/// </ul>
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
/// <li> Warning step size: 1000.
/// <li> Show period: 25. 
/// </ul>
///
/// @see OptimizationAlgorithm.

ConjugateGradient::ConjugateGradient(void) : OptimizationAlgorithm()
{
   // Optimization operators

   searchDirectionMethod = PolakRibiere;
   optimalStepSizeMethod = BrentMethod;

   // Optimization parameters

   firstStepSize = 1.0e-3;
   optimalStepSizeTolerance = 1.0e-3;

   // Stopping criteria

   evaluationGoal = -1.0e69;
   gradientNormGoal = 0.0;
   maximumTime = 1.0e6;
   maximumNumberOfIterations = 1000;

   // User stuff

   warningStepSize = 1000.0;
   showPeriod = 25;
}


// DESTRUCTOR

/// Destructor.

ConjugateGradient::~ConjugateGradient(void)
{

}


// METHODS

// SearchDirectionMethod getSearchDirectionMethod(void) method

/// This method returns the search direction method used for optimization.
///
/// @see getFletcherReevesSearchDirection(Vector<double>, Vector<double>, Vector<double>).
/// @see getPolakRibiereSearchDirection(Vector<double>, Vector<double>, Vector<double>).

ConjugateGradient::SearchDirectionMethod 
ConjugateGradient::getSearchDirectionMethod(void)
{
   return(searchDirectionMethod);
}


// OptimalStepSizeMethod getOptimalStepSizeMethod(void) method

/// This method returns the optimal step size method used for optimization.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

ConjugateGradient::OptimalStepSizeMethod
ConjugateGradient::getOptimalStepSizeMethod(void)
{
   return(optimalStepSizeMethod);
}


// Vector<double> getInitialArgument(void)

/// This method returns the initial objective function argument to be used by 
/// the conjugate gradient method for optimization. 

Vector<double> ConjugateGradient::getInitialArgument(void)
{
   return(initialArgument);
}


// double getGradientNormGoal(void) method

/// This method returns the goal value for the norm of the objective function
/// gradient.
/// This is used as a stopping criterium when optimizing a function.

double ConjugateGradient::getGradientNormGoal(void)
{
   return(gradientNormGoal);
}


// int getMaximumNumberOfIterations(void) method

/// This method returns the maximum number of iterations to be performed by the 
/// conjugate gradient method during the optimization process. 
/// This is used as a stopping criterium when optimizing an objective function.

int ConjugateGradient::getMaximumNumberOfIterations(void)
{
   return(maximumNumberOfIterations);
}


// double getFirstStepSize(void) method

/// This method returns the initial step size for line search
/// in the first iteration of the conjugate gradient.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

double ConjugateGradient::getFirstStepSize(void)
{
   return(firstStepSize);
}


// double getOptimalStepSizeTolerance(void) method

/// This method returns the tolerance value in line search.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

double ConjugateGradient::getOptimalStepSizeTolerance(void)
{
   return(optimalStepSizeTolerance);
}


// double getWarningStepSize(void) method

/// This method returns the step size value at wich a warning message is
/// written to the screen during line search.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

double ConjugateGradient::getWarningStepSize(void)
{
   return(warningStepSize);
}


// int getShowPeriod(void) method

/// This method returns the number of iterations between the optimization 
/// showing progress. 

int ConjugateGradient::getShowPeriod(void)
{
   return(showPeriod);    
}


// void setSearchDirectionMethod(SearchDirectionMethod) method

/// This method sets a new search direction method to be used for otpimization
/// with the conjugate gradient method.
///
/// @param newSearchDirectionMethod: Search direction method.
///
/// @see getFletcherReevesSearchDirection(Vector<double>, Vector<double>, Vector<double>).
/// @see getPolakRibiereSearchDirection(Vector<double>, Vector<double>, Vector<double>).

void ConjugateGradient::setSearchDirectionMethod
(ConjugateGradient::SearchDirectionMethod newSearchDirectionMethod)
{
   searchDirectionMethod = newSearchDirectionMethod;
}


// void setOptimalStepSizeMethod(StepSizeMethod) method

/// This method sets a new optimal step size method to be used for optimization
/// with the conjugate gradient method.
///
/// @param newOptimalStepSizeMethod: Optimal step size method.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

void ConjugateGradient::setOptimalStepSizeMethod
(ConjugateGradient::OptimalStepSizeMethod newOptimalStepSizeMethod)
{
   optimalStepSizeMethod = newOptimalStepSizeMethod;
}


// void setInitialArgument(Vector<double>) method

/// This method sets a new initial objective function argument to be used by 
/// the conjugate gradient method for optimization. 
///
/// @param newInitialArgument: Initial argument Vector.

void ConjugateGradient::setInitialArgument(Vector<double> newInitialArgument)
{
   int size = newInitialArgument.getSize();

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   if(size != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class. "
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

void ConjugateGradient::setGradientNormGoal(double newGradientNormGoal)
{
   if(gradientNormGoal < 0.0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class." << std::endl
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

/// This method sets a new maximum number of conjugate gradient iterations for
/// optimization.  
///
/// @param newMaximumNumberOfIterations: Maximum number of iterations.

void ConjugateGradient
::setMaximumNumberOfIterations(int newMaximumNumberOfIterations)
{
   if(newMaximumNumberOfIterations <= 0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class." << std::endl
                << "void setMaximumNumberOfIterations(int) method."
                << std::endl
                << "Maximum number of iterations must be greater than 0."
                << std::endl
                << std::endl;

      exit(1);
   }

   maximumNumberOfIterations = newMaximumNumberOfIterations;
}


// void setShowPeriod(int) method

/// This sets a new number of iterations between the optimization showing progress. 
///
/// @param newShowPeriod: Show period.

void ConjugateGradient::setShowPeriod(int newShowPeriod)
{
   if(newShowPeriod <= 0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class." << std::endl
                << "void setShowPeriod(int) method."
                << std::endl
                << "Show period must be greater than 0."
                << std::endl << std::endl;

      exit(1);
   }
   else
   {
      showPeriod = newShowPeriod;
   }
}


// void setFirstStepSize(double) method

/// This method sets a new initial step size for line search
/// in the first iteration of the conjugate gradient.
///
/// @param newFirstStepSize: First step size.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

void ConjugateGradient::setFirstStepSize(double newFirstStepSize)
{ 
   if(newFirstStepSize <= 0.0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class." << std::endl
                << "void setFirstStepSize(double) method."
                << std::endl
                << "First step size must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   // Set first step size

   firstStepSize = newFirstStepSize;
}


// void setToleranceInOptimalStepSize(double) method

/// This method sets a new tolerance value to be used in line line search.
///
/// @param newOptimalStepSizeTolerance: Tolerance value.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

void ConjugateGradient
::setOptimalStepSizeTolerance(double newOptimalStepSizeTolerance)
{
   if(optimalStepSizeTolerance <= 0.0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class. "
                << "void setOptimalStepSizeTolerance(double) method."
                << std::endl
                << "Tolerance must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      optimalStepSizeTolerance = newOptimalStepSizeTolerance;
   }
}


// void setWarningStepSize(double) method

/// This method sets a new step size value at wich a warning message is written
/// to the screen during line search.
///
/// @param newWarningStepSize: Warning step size value.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).

void ConjugateGradient::setWarningStepSize(double newWarningStepSize)
{
   if(newWarningStepSize <= 0.0)
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class. " << std::endl
                << "void setWarningStepSize(double) method." << std::endl
                << "Warning step size must be greater than 0." << std::endl
                << std::endl;

      exit(1);
   }

   warningStepSize = newWarningStepSize;
}


// double getFletcherReevesParameter(Vector<double>, Vector<double>) method

/// This method returns the Fletcher-Reeves parameter used to get the search
/// direction.
///
/// @param oldGradient: Previous objective function gradient.
/// @param gradient: Current objective function gradient.
///
/// @see getFletcherReevesSearchDirection(Vector<double>, Vector<double>, Vector<double>).

double ConjugateGradient
::getFletcherReevesParameter(Vector<double> oldGradient, Vector<double> gradient)
{
   double FletcherReevesParameter = 0.0;

   int n = gradient.getSize();

   double numerator = 0.0;
   double denominator = 0.0;

   for (int i = 0; i < n; i++)
   {
      numerator += pow(gradient[i], 2);

      denominator += pow(oldGradient[i], 2);
   }

   // Prevent a possible division by 0

   if (denominator == 0.0)
   {
      FletcherReevesParameter = 0.0;
   }
   else
   {
      FletcherReevesParameter = numerator/denominator;
   }

   // Bound the Fletcher-Reeves parameter between 0 and 1

   if (FletcherReevesParameter < 0.0)
   {
      FletcherReevesParameter = 0.0;
   }
   else if (FletcherReevesParameter > 1.0)
   {
      FletcherReevesParameter = 1.0;
   }

   return(FletcherReevesParameter);
}


// double getPolakRibiereParameter(Vector<double>, Vector<double>) method

/// This method returns the Polak-Ribiere parameter used to get the search
/// direction.
///
/// @param oldGradient: Previous objective function gradient.
/// @param gradient: Current objective function gradient.
///
/// @see getPolakRibiereSearchDirection(Vector<double>, Vector<double>, Vector<double>).

double ConjugateGradient
::getPolakRibiereParameter(Vector<double> oldGradient, Vector<double> gradient)
{
   int n = gradient.getSize();

   double PolakRibiereParameter = 0.0;

   double numerator = 0.0;
   double denominator = 0.0;

   for (int i = 0; i < n; i++)
   {
      numerator += (gradient[i] - oldGradient[i])*gradient[i];
      denominator += pow(oldGradient[i], 2);
   }

   // Prevent a possible division by 0

   if (denominator == 0.0)
   {
      PolakRibiereParameter = 0.0;
   }
   else
   {
      PolakRibiereParameter = numerator/denominator;
   }

   // Bound the Polak-Ribiere parameter between 0 and 1

   if (PolakRibiereParameter < 0.0)
   {
      PolakRibiereParameter = 0.0;
   } 
   else if (PolakRibiereParameter > 1.0)
   {
      PolakRibiereParameter = 1.0;
   }

   return(PolakRibiereParameter);
}


// Vector<double> getFletcherReevesSearchDirection
// (Vector<double>, Vector<double>, Vector<double>) method

/// This method returns the search direction using the Fletcher-Reeves
/// update.
///
/// @param oldGradient: Previous objective function gradient.
/// @param gradient: Current objective function gradient.
/// @param oldSearchDirection: Previous search direction Vector.
///
/// @see getFletcherReevesParameter(Vector<double>, Vector<double>).
/// @see getPolakRibiereParameter(Vector<double>, Vector<double>)
/// @see getPolakRibiereSearchDirection(Vector<double>, Vector<double>, Vector<double>)

Vector<double> ConjugateGradient::getFletcherReevesSearchDirection
(Vector<double> oldGradient, Vector<double> gradient, Vector<double> oldSearchDirection)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> FletcherReevesSearchDirection(numberOfVariables, 0.0);

   double FletcherReevesParameter = getFletcherReevesParameter(oldGradient, gradient);

   for (int i = 0; i < numberOfVariables; i++)
   {
      FletcherReevesSearchDirection[i]
      = -1.0*gradient[i] + FletcherReevesParameter*oldSearchDirection[i];
   }

   return(FletcherReevesSearchDirection);
}


// Vector<double> getPolakRibiereSearchDirection(Vector<double>, Vector<double>,
// Vector<double>) method

/// This method returns the search direction using the Polak-Ribiere
/// update.
///
/// @param oldGradient: Previous objective function gradient.
/// @param gradient: Current objective function gradient.
/// @param oldSearchDirection: Previous search direction Vector.
///
/// @see getPolakRibiereParameter(Vector<double>, Vector<double>).
/// @see getFletcherReevesParameter(Vector<double>, Vector<double>).
/// @see getFletcherReevesSearchDirection(Vector<double>, Vector<double>, Vector<double>).

Vector<double> ConjugateGradient::getPolakRibiereSearchDirection
(Vector<double> oldGradient, Vector<double> gradient, Vector<double> oldSearchDirection)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> PolakRibiereSearchDirection(numberOfVariables, 0.0);

   double PolakRibiereParameter = getPolakRibiereParameter(oldGradient, gradient);

   for (int i = 0; i < numberOfVariables; i++)
   {
      PolakRibiereSearchDirection[i]
      = -1.0*gradient[i] + PolakRibiereParameter*oldSearchDirection[i];
   }

   return(PolakRibiereSearchDirection);
}


// double getGoldenSectionOptimalStepSize(double, Vector<double>, Vector<double>, double)
// method

/// This method returns the optimal step size by searching in a given
/// direction to locate the minimum of the objective function in that
/// direction. It uses the golden section method.
///
/// @param initialStepSize: Initial step size in line search.
/// @param evaluation: Objective function evaluation value.
/// @param argument: Objective function argument Vector.
/// @param searchDirection: Search direction Vector.
///
/// @see getBrentMethodOptimalStepSize(double, double, Vector<double>, Vector<double>).
 
double ConjugateGradient::getGoldenSectionOptimalStepSize
(double initialStepSize, double evaluation,
Vector<double> argument, Vector<double> searchDirection)
{
   double optimalStepSize = 0.0;

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> potentialArgument(numberOfVariables, 0.0);

   double a = 0.0;
   double evaluationA = 0.0;
   double b = 0.0;
   double evaluationB = 0.0;

   double c = 0.0;	
   double evaluationC = 0.0;
   double d = 0.0;	
   double evaluationD = 0.0;

   double tau = (3.0-sqrt(5.0))/2.0; // 0.382

   // Start golden section search

   // Set initial a point

   a = 0.0;

   // Calculate evaluation for a

   evaluationA = evaluation;

   // Set initial b point

   b = initialStepSize;

   // Calculate evaluation for b

   for (int i = 0; i < numberOfVariables; i++)
   {
      potentialArgument[i] = argument[i] + searchDirection[i]*b;
   }

   evaluationB = objectiveFunction->getEvaluation(potentialArgument);

   // Find initial interval where minimum evaluation occurs

   while(evaluationA > evaluationB)
   {
      // Set new b

      b = 2.0*b;

      if(b >= warningStepSize)
      {
         std::cout << std::endl
                   << "Warning: Step size is " << b
                   << std::endl;
      }

      // Calculate evaluation for new b

      for (int i = 0; i < numberOfVariables; i++)
      {
         potentialArgument[i] = argument[i] + searchDirection[i]*b;
      }

      evaluationB = objectiveFunction->getEvaluation(potentialArgument);
   }

   // Initialize c and d (interior points for line minimization)

   // Initialize c point

   c = a + tau*(b-a);

   // Calculate evaluation for c

   for (int i = 0; i < numberOfVariables; i++)
   {
      potentialArgument[i] = argument[i] + searchDirection[i]*c;
   }

   evaluationC = objectiveFunction->getEvaluation(potentialArgument);

   // Initialize d point

   d = b - tau*(b-a);

   // Calculate evaluation for d

   for (int i = 0; i < numberOfVariables; i++)
   {
      potentialArgument[i] = argument[i] + searchDirection[i]*d;
   }

   evaluationD = objectiveFunction->getEvaluation(potentialArgument);

   // Reduce the interval with the golden section algorithm

   while(b-a > optimalStepSizeTolerance)
   {
      Vector<double> evaluationVectorLeft(3, 0.0);
      evaluationVectorLeft[0] = evaluationA;
      evaluationVectorLeft[1] = evaluationC;
      evaluationVectorLeft[2] = evaluationD;

      double minimumEvaluationLeft = getMinimum(evaluationVectorLeft);

      Vector<double> evaluationVectorRight(3, 0.0);
      evaluationVectorRight[0] = evaluationB;
      evaluationVectorRight[1] = evaluationC;
      evaluationVectorRight[2] = evaluationD;

      double minimumEvaluationRight = getMinimum(evaluationVectorRight);

      if((evaluationC <= evaluationD && evaluationB >= minimumEvaluationLeft)
      || (evaluationA <=  minimumEvaluationRight))

      // There is a minimum between a and b
      {
         b=d; 
         d=c; 

         evaluationB = evaluationD;
         evaluationD = evaluationC;

         // Set new c point

         c = a + tau*(b-a);

         // Calculate evaluation for new c

         for (int i = 0; i < numberOfVariables; i++)
         {
            potentialArgument[i] = argument[i] + searchDirection[i]*c;
         }

         evaluationC = objectiveFunction->getEvaluation(potentialArgument);
      }
      else if((evaluationD <= evaluationC && evaluationA >= minimumEvaluationRight)
      || (evaluationB <= minimumEvaluationLeft))

      // There is a minimum between c and b
      {
         a = c; 
         c = d; 

         evaluationA = evaluationC;
         evaluationC = evaluationD;

         // Set new d point

         d = b - tau*(b-a);

         // Calculate evaluation for new d

         for (int i = 0; i < numberOfVariables; i++)
         {
            potentialArgument[i] = argument[i] + searchDirection[i]*d;
         }

         evaluationD = objectiveFunction->getEvaluation(potentialArgument);
      }
      else
      {
         std::cout << std::endl
                   << "Error: ConjugateGradient class." << std::endl
                   << "double getGoldenSectionStepSize "
                   << "(double, double, Vector<double>, Vector<double>, double) method."
                   << std::endl
                   << "Unable to find were the minimum is." << std::endl;

         exit(1);
      }
   }

   // Get minimum evaluation and optimal step size among points A, B, C and D

   double minimumEvaluation = evaluation;

   if(evaluationA < minimumEvaluation)
   {
      minimumEvaluation  = evaluationA;
      optimalStepSize = a;
   }
   else if(evaluationB < minimumEvaluation)
   {
      minimumEvaluation = evaluationB;
      optimalStepSize = b;
   }
   else if(evaluationC < minimumEvaluation)
   {
      minimumEvaluation = evaluationC;
      optimalStepSize = c;
   }
   else if(evaluationD < minimumEvaluation)
   {
      minimumEvaluation = evaluationD;
      optimalStepSize = d;
   }

   return(optimalStepSize);
}


// double getBrentMethodOptimalStepSize(double, Vector<double>, Vector<double>, double)
// method

/// This method returns the optimal step size by searching in a given
/// direction to locate the minimum of the objective function in that
/// direction. It uses the Brent's method.
///
/// @param initialStepSize: Initial step size in line search.
/// @param evaluation: Objective function evaluation value.
/// @param argument: Objective function argument Vector.
/// @param searchDirection: Search direction Vector.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>).

double ConjugateGradient::getBrentMethodOptimalStepSize
(double initialStepSize, double evaluation,
Vector<double> argument, Vector<double> searchDirection)
{
   double optimalStepSize = 0.0;

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> potentialArgument(numberOfVariables, 0.0);

   double a = 0.0;
   double evaluationA = 0.0;
   double b = 0.0;
   double evaluationB = 0.0;

   double u = 0.0;
   double evaluationU = 0.0;
   double v = 0.0;
   double evaluationV = 0.0;
   double w = 0.0;
   double evaluationW = 0.0;
   double x = 0.0;
   double evaluationX = 0.0;

   double tau = (3.0-sqrt(5.0))/2.0; // 0.382

   // Set initial a point

   a = 0.0;

   // Calculate evaluation for a

   evaluationA = evaluation;

   // Set initial b point

   b = initialStepSize;

   // Calculate evaluation for b

   for (int i = 0; i < numberOfVariables; i++)
   {
      potentialArgument[i] = argument[i] + searchDirection[i]*b;
   }

   evaluationB = objectiveFunction->getEvaluation(potentialArgument);

   // Find initial interval where minimum evaluation occurs

   while(evaluationA > evaluationB)
   {
      // Set new b

      b = 2.0*b;

      if(b >= warningStepSize)
      {
         std::cout << std::endl
                   << "Warning: Step size is  " << b
                   << std::endl;
      }

      // Calculate evaluation for new b

      for (int i = 0; i < numberOfVariables; i++)
      {
         potentialArgument[i] = argument[i] + searchDirection[i]*b;
      }

      evaluationB = objectiveFunction->getEvaluation(potentialArgument);
   }

   // Get intermediate point V

   v = a + tau*(b-a);

   // Calculate evaluation for V

   for (int i = 0; i < numberOfVariables; i++)
   {
      potentialArgument[i] = argument[i] + searchDirection[i]*v;
   }

   evaluationV = objectiveFunction->getEvaluation(potentialArgument);

   // Set initial W and X points

   w = v;
   evaluationW = evaluationV;

   x = v;
   evaluationX = evaluationV;

   // Maximum and minimum intervals ???

   double intervalLength = 0.0;

   bool goldenSection = false;

   // Reduce the interval

   while(b-a > optimalStepSizeTolerance)
   {
      // Quadratic interpolation

      if(w != x && w != v && x != v) // Can construct parabola
      {
         // zz vector

         Vector<double> stepSizeVector(3, 0.0);
         stepSizeVector[0] = v;
         stepSizeVector[1] = w;
         stepSizeVector[2] = x;

         std::sort(stepSizeVector.begin(), stepSizeVector.end(), std::less<double>());

         // pp vector

         Vector<double> evaluationVector(3, 0.0);

         for(int i = 0; i < 3; i++)
         {
            if(stepSizeVector[i] == v)
            {
               evaluationVector[i] = evaluationV;
            }
            else if(stepSizeVector[i] == w)
            {
               stepSizeVector[i] = evaluationW;
            }
            else if(stepSizeVector[i] == x)
            {
               stepSizeVector[i] = evaluationX;
            }
            else
            {
               std::cout << std::endl
                         << "Error: ConjugateGradient class." << std::endl
                         << "double getBrentMethodOptimalStepSize" << std::endl
                         << "(double, double, Vector<double>, Vector<double>) method."
                         << std::endl
                         << "Unable to construct step size and evaluation vectors right."
                         << std::endl << std::endl;

               exit(1);
            }
         }

         // xStar is the minimum of the parabola through the three step size points

         double numerator
         = (pow(stepSizeVector[2],2) - pow(stepSizeVector[1],2))*evaluationVector[0]
         + (pow(stepSizeVector[1],2) - pow(stepSizeVector[0],2))*evaluationVector[2]
         + (pow(stepSizeVector[0],2) - pow(stepSizeVector[2],2))*evaluationVector[1];

         double denominator
         = (stepSizeVector[2] - stepSizeVector[1])*evaluationVector[0]
         + (stepSizeVector[1] - stepSizeVector[0])*evaluationVector[2]
         + (stepSizeVector[0] - stepSizeVector[2])*evaluationVector[1];

         double xStar = 0.5*numerator/denominator;

         if(xStar < b && a < xStar) // xStar is in [a,b]
         {
            u = xStar;

            // Good, no need to perform golden section

            goldenSection = false;
         }
         else // xStar is not in [a,b]
         {
            // Bad, need to perform golden section

            goldenSection = true;
         }
      }
      else // Cannot construct parabola
      {
         // Bad, need to perform golden section

         goldenSection = true;
      }

      //
      // Golden section
      //

      if(goldenSection == true)
      {
         if(x > (a+b)/2.0)
         {
            u = x-tau*(x-a);
         }
         else
         {
            u = x+tau*(b-x);
         }
      }

      // Calculate evaluation for U

      for (int i = 0; i < numberOfVariables; i++)
      {
         potentialArgument[i] = argument[i] + searchDirection[i]*u;
      }

      evaluationU = objectiveFunction->getEvaluation(potentialArgument);

      // Update points

      if(evaluationU < evaluationX)
      {
         if(u < x)
         {
            b = x;
            evaluationB = evaluationX;
         }
         else
         {
            a = x;
            evaluationA = evaluationX;
         }

         v = w;
         evaluationV = evaluationW;

         w = x;
         evaluationW = evaluationX;

         x = u;
         evaluationX = evaluationU;
      }
      else
      {
         if(u < x)
         {
            a = u;
            evaluationA = evaluationU;
         }
         else
         {
            b = u;
            evaluationB = evaluationU;
         }

         if((evaluationU < evaluationW) || (w == x))
         {
             v = w;
             evaluationV = evaluationW;

             w = u;
             evaluationW = evaluationU;
         }
         else if((evaluationU < evaluationV) || (v == x) || (v == w))
         {
            v = u;
            evaluationV = evaluationU;
         }
      }
   } // while loop

   // Get minimum evaluation and optimal step size among points A, B, V, W and X

   double minimum = evaluation;

   if(evaluationA < minimum)
   {
      minimum = evaluationA;
      optimalStepSize = a;
   }
   else if(evaluationB < minimum)
   {
      minimum = evaluationB;
      optimalStepSize = b;
   }
   else if(evaluationV < minimum)
   {
      minimum = evaluationV;
      optimalStepSize = v;
   }
   else if(evaluationW < minimum)
   {
      minimum = evaluationW;
      optimalStepSize = w;
   }
   else if(evaluationX < minimum)
   {
      minimum = evaluationX;
      optimalStepSize = x;
   }

   return(optimalStepSize);
}


// void getMinimalArgument(void) method

/// This method optimizes an objective function according to the 
/// conjugate gradient algorithm. 
/// It returns the minimal argument of the objective function.
/// Optimization occurs according to the optimization operators and the 
/// optimization parameters.

Vector<double> ConjugateGradient::getMinimalArgument(void)
{
   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   Vector<double> minimalArgument(numberOfVariables, 0.0);
   Vector<double> argument(numberOfVariables, 0.0);

   // Evaluation history vector

   Vector<double> newEvaluationHistory(maximumNumberOfIterations+1, 0.0);
   evaluationHistory = newEvaluationHistory;

   // Objective function gradient norm optimization history vector

   Vector<double> newGradientNormHistory(maximumNumberOfIterations+1, 0.0);
   gradientNormHistory = newGradientNormHistory;

   double evaluation = 0.0;

   Vector<double> gradient(numberOfVariables, 0.0);
   Vector<double> oldGradient(numberOfVariables, 0.0);

   double gradientNorm = 0.0;

   Vector<double> searchDirection(numberOfVariables, 0.0);
   Vector<double> oldSearchDirection(numberOfVariables, 0.0);

   double slope = 0.0;

   double initialStepSize = 0.0;
   double optimalStepSize = 0.0;
   double oldOptimalStepSize = 0.0;

   time_t beginningTime, currentTime;
   double elapsedTime = 0.0;

   // Set beginning optimization time 

   time(&beginningTime);

   std::cout << std::endl
             << "Getting minimal argument with conjugate gradient..." 
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

      objectiveFunction->print(argument);
            
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

      objectiveFunction->print(argument);
              
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

   objectiveFunction->print(argument);

   // Loop over iterations

   for(int iteration = 1; iteration <= maximumNumberOfIterations; iteration++)
   {

      // Get search direction

      if(iteration % numberOfVariables == 0)
      {
         // Search direction is set to negative gradient

         for(int i = 0; i < numberOfVariables; i++)
         {
            searchDirection[i] = -1.0*gradient[i];
         } 
      }
      else
      {
         // Conjugate gradient search direction

         switch(searchDirectionMethod)
         {
            case FletcherReeves:

               searchDirection = getFletcherReevesSearchDirection
               (oldGradient, gradient, oldSearchDirection);

               break;

            case PolakRibiere:

              searchDirection = getPolakRibiereSearchDirection
              (oldGradient, gradient, oldSearchDirection);

               break;
         }
      }

      // Calculate slope

      slope = 0.0;

      for(int i = 0; i < numberOfVariables; i++)
      {
         slope += gradient[i]*searchDirection[i];
      }

      // Check for a descent direction 

      if(slope >= 0.0)
      {
         std::cout << std::endl
                   << "Warning: Function slope is equal or greater than zero."
                   << std::endl
                   << "Search direction reset to negative gradient."
                   << std::endl;

         // Reset search direction

         for(int i = 0; i < numberOfVariables; i++)
         {
            searchDirection[i] = -1.0*gradient[i];
         }

         // Calculate new slope

         slope = 0.0;

         for(int i = 0; i < numberOfVariables; i++)
         {
            slope += gradient[i]*searchDirection[i];
         }
      }

      // Get optimal step size

      if(iteration == 1)
      {
         initialStepSize = firstStepSize;
      }
      else
      {
         initialStepSize = oldOptimalStepSize;
      }

      switch(optimalStepSizeMethod)
      {
         case GoldenSection:

            optimalStepSize = getGoldenSectionOptimalStepSize
            (initialStepSize, evaluation, argument, searchDirection);


            break;

         case BrentMethod:

            optimalStepSize = getBrentMethodOptimalStepSize
            (initialStepSize, evaluation, argument, searchDirection);

            break;
      }

      // Optimal step size stopping criterium

      if(optimalStepSize == 0.0)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Optimal step size is zero." << std::endl;

         std::cout << "Final evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Final gradient norm: " << gradientNorm << ";" << std::endl;

         objectiveFunction->print(argument);

         break;
      }

      // Get new free parameters

      for (int i = 0; i < numberOfVariables; i++)
      {
         argument[i] += optimalStepSize*searchDirection[i];
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

            std::cout << "Evaluation: " << evaluation << ";" << std::endl;
            std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;

            objectiveFunction->print(argument);

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

         objectiveFunction->print(argument);

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

         objectiveFunction->print(argument);

         break;
      }

      // Maximum number of iterations

      if (iteration == maximumNumberOfIterations)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << ": "
                   << "Maximum number of iterations reached."
                   << std::endl;

         std::cout << "Evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;

         objectiveFunction->print(argument);

         break;
      }

      // Progress

      if(iteration % showPeriod == 0)
      {
         std::cout << std::endl
                   << "Iteration " << iteration << "; " << std::endl;

         std::cout << "Evaluation: " << evaluation << ";" << std::endl;
         std::cout << "Gradient norm: " << gradientNorm << ";" << std::endl;
         
         objectiveFunction->print(argument);
      }

      // Update

      oldGradient = gradient;
      oldSearchDirection = searchDirection;
      oldOptimalStepSize = optimalStepSize;
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

/// This method prints to the screen the optimization operators and
/// the optimization parameters concerning the conjugate gradient object:
///
/// Initial argument.
///
/// Optimization operators:
/// <ul>
/// <li> Search direction method.
/// <li> Optimal step size method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> First step size.
/// <li> Optimal step size tolerance.
/// </ul>
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
/// <li> Warning step size.
/// <li> Show period. 
/// </ul>

void ConjugateGradient::print(void)
{
   std::cout << std::endl
             << "Conjugate Gradient Object." << std::endl; 

   int numberOfVariables = objectiveFunction->getNumberOfVariables();

   // Initial argument

   std::cout << std::endl
             << "Initial argument:" << std::endl;

   for(int i = 0; i < numberOfVariables; i++)
   {
      std::cout << initialArgument[i] << " ";        
   }
   
   std::cout << std::endl;

   // Optimization operators

   std::cout << std::endl
             << "Optimization operators:" << std::endl;

   // Search direction method

   std::cout << "Search direction method:" << std::endl;

   switch(searchDirectionMethod)
   {
      case FletcherReeves:

      std::cout << "Fletcher-Reeves" << std::endl;

      break;

      case PolakRibiere:

      std::cout << "Polak-Ribiere" << std::endl;

      break;
   }


   // Optimal step size method

   std::cout << "Optimal step size method:" << std::endl;

   switch(optimalStepSizeMethod)
   {
      case GoldenSection:

         std::cout << "Golden section" << std::endl;

      break;

      case BrentMethod:

         std::cout << "Brent Method" << std::endl;

      break;
   }


   // Optimization parameters

   std::cout << std::endl
             << "Optimization parameters: " << std::endl
             << "First step size: " << std::endl
             << firstStepSize << std::endl
             << "Optimal step size tolerance: " << std::endl
             << optimalStepSizeTolerance << std::endl;

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
             << "Warning step size: " << std::endl
             << warningStepSize << std::endl
             << "Show period: " << std::endl
             << showPeriod
             << std::endl;
}


// void save(char*) method

/// This method saves the conjugate gradient object to a data file. 
///
/// Initial argument.
///
/// Optimization operators:
/// <ul>
/// <li> Search direction method.
/// <li> Optimal step size method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> First step size.
/// <li> Optimal step size tolerance.
/// </ul>
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
/// <li> Warning step size.
/// <li> Show period. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see load(char*).

void ConjugateGradient::save(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: ConjugateGradient class. " << std::endl
                << "void save(char*) method." << std::endl
                << "Cannot open conjugate gradient object data file."  << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Saving conjugate gradient object to data file..." << std::endl;
   }

   // Write file header

   file << "% Purple: An Open Source Numerical Optimization C++ Library." 
        << std::endl 
        << "% Conjugate Gradient Object." << std::endl; 

   int numberOfVariables = objectiveFunction->getNumberOfVariables();
 
   // Initial argument

   file << "InitialArgument:" << std::endl;

   for(int i = 0; i < numberOfVariables; i++)
   {
      file << initialArgument[i] << " ";        
   }
   
   file << std::endl;

   // Optimization operators

   // Search direction method

   file << "SearchDirectionMethod:" << std::endl;


   switch(searchDirectionMethod)
   {
      case FletcherReeves:

      file << "FletcherReeves" << std::endl;

      break;

      case PolakRibiere:

      file << "PolakRibiere" << std::endl;

      break;
   }


   // Optimal step size method

   file << "OptimalStepSizeMethod:" << std::endl;

   switch(optimalStepSizeMethod)
   {
      case GoldenSection:

      file << "GoldenSection" << std::endl;

      break;

      case BrentMethod:

      file << "BrentMethod" << std::endl;

      break;
   }


   // Optimization parameters

   file << "FirstStepSize:" << std::endl
        << firstStepSize << std::endl
        << "OptimalStepSizeTolerance:" << std::endl
        << optimalStepSizeTolerance << std::endl;

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

   file << "WarningStepSize: " << std::endl
        << warningStepSize << std::endl
        << "ShowPeriod: " << std::endl
        << showPeriod << std::endl;

   file.close();
}


// void load(char*) method

/// This method loads a conjugate gradient object from a data file. 
/// Please mind about the file format, wich is specified in the User's Guide. 
///
///
/// Initial argument.
///
/// Optimization operators:
/// <ul>
/// <li> Search direction method.
/// <li> Optimal step size method.
/// </ul>
///
/// Optimization parameters:
/// <ul>
/// <li> First step size.
/// <li> Optimal step size tolerance.
/// </ul>
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
/// <li> Warning step size.
/// <li> Show period. 
/// </ul>
///
/// @param filename: Filename.
///
/// @see save(char*).

void ConjugateGradient::load(char* filename)
{
   // File

   std::fstream file;

   file.open(filename, std::ios::in);

   if(!file.is_open())
   {
      std::cout << std::endl
                << "Error: ConjugateGradient class." << std::endl
                << "void load(char*) method." << std::endl
                << "Cannot open conjugate gradient object data file."  << std::endl;

      exit(1);
   }
   else
   {
      std::cout << std::endl
                << "Loading conjugate gradient object from data file..."  << std::endl;
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

   // Optimization operators

   // Search direction method

   file >> word;

   if(word == "FletcherReeves")
   {
      searchDirectionMethod = FletcherReeves;
   }
   else if(word == "PolakRibiere")
   {
      searchDirectionMethod = PolakRibiere;
   }

   // Optimal step size method

   file >> word;

   file >> word;

   if(word == "GoldenSection")
   {
      optimalStepSizeMethod = GoldenSection;
   }
   else if(word == "BrentMethod")
   {
      optimalStepSizeMethod = BrentMethod;
   }

   // Optimization parameters

   // First step size

   file >> word;

   file >> firstStepSize;

   // Tolerance

   file >> word;

   file >> optimalStepSizeTolerance;

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

   // Warning step size

   file >> word;

   file >> warningStepSize;

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

void ConjugateGradient::saveOptimizationHistory(char* filename)
{
   std::fstream file;

   file.open(filename, std::ios::out);

   if(!file.is_open())
   {
      std::cout << std::endl 
                << "Error: ConjugateGradient class. "
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
        << "% Conjugate Gradient Optimization History." << std::endl
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


// double getMinimum(Vector<double>) method

/// This method returns the minimum element in a Vector of double precision numbers.
///
/// @param v: Vector of double precision numbers.
///
/// @see getGoldenSectionStepSize(double, double, Vector<double>, Vector<double>)

double ConjugateGradient::getMinimum(Vector<double> v)
{
   int n = v.getSize();

   double minimum = v[0];

   for(int i = 1; i < n; i++)
   {
      if(v[i] < minimum)
      {
         v[i] = minimum;
      }
   }

   return(minimum);
}


// double getMaximum(Vector<double>) method

/// This method returns the maximum element in a Vector of double precision numbers.
///
/// @param v: Vector of double precision numbers.
///
/// @see getGoldenSectionOptimalStepSize(double, double, Vector<double>, Vector<double>)

double ConjugateGradient::getMaximum(Vector<double> v)
{
   int n = v.getSize();

   double maximum = v[0];

   for(int i = 1; i < n; i++)
   {
      if(v[i] > maximum)
      {
         v[i] = maximum;
      }
   }

   return(maximum);
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
