/******************************************************************************/
/*                                                                            */
/*   O B J E C T I V E   F U N C T I O N   C L A S S                          */
/*                                                                            */
/*   Roberto Lopez                                                            */
/*   International Center for Numerical Methods in Engineering (CIMNE)        */
/*   Technical University of Catalonia (UPC)                                  */
/*   Barcelona, Spain                                                         */
/*   E-mail: rlopez@cimne.upc.edu                                             */
/*                                                                            */
/******************************************************************************/


#include "ObjectiveFunction.h"

#include<iostream>
#include <math.h>
#include <sstream>
#include <stdexcept>

namespace Purple
{

// GENERAL CONSTRUCTOR

/// General constructor. It creates an objective function object. 
/// It also initializes all the rest of class members to their default values:
///
/// <ul>
/// <li> Number of evaluations = 0.
/// <li> Epsilon: 1.0e-6.
/// </ul> 

ObjectiveFunction::ObjectiveFunction(void)
{
   numberOfEvaluations = 0;

   epsilon = 1.0e-6;
}


// DESTRUCTOR

/// Destructor.

ObjectiveFunction::~ObjectiveFunction(void)
{

}


// METHODS

// int getNumberOfVariables(void) method

/// This method returns the number of variables in the objective function.

int ObjectiveFunction::getNumberOfVariables(void)
{
   return(numberOfVariables);
}


// Vector<double> getLowerBound(void)

/// This method returns the lower bound of the objective function domain.
///
/// @see getUpperBound(void)
/// @see getDomain(void)

Vector<double> ObjectiveFunction::getLowerBound(void)
{
   return(lowerBound);               
}

// Vector<double> getUpperBound(void)

/// This method returns the upper bound of the objective function domain.
///
/// @see getLowerBound(void)
/// @see getDomain(void)

Vector<double> ObjectiveFunction::getUpperBound(void)
{
   return(upperBound);               
}

// Matrix<double> getDomain(void)

/// This method returns the objective function domain in a single Matrix.
/// Row 0 contains the lower bound of the domain.
/// Row 1 contains the upper bound of the domain.
///
/// @see getLowerBound(void)
/// @see getUpperBound(void)

Matrix<double> ObjectiveFunction::getDomain(void)
{
   Matrix<double> domain(2, numberOfVariables, 0.0);
   
   for(int i = 0; i< numberOfVariables; i++)
   {
      domain[0][i] = lowerBound[i];
      domain[1][i] = upperBound[i];              
   }

   return(domain);               
}


// double getEpsilon(void) method

/// This method returns the epsilon value to be used for  numerical 
/// differentiation.
///
/// @see getGradient(void)
/// @see getHessian(void)

double ObjectiveFunction::getEpsilon(void)
{
   return(epsilon);
}


// int getNumberOfEvaluations(void) method

/// This method returns the number of calls to the getEvaluation(Vector<double>) 
/// method.
///
/// @see getEvaluation(Vector<double>).

int ObjectiveFunction::getNumberOfEvaluations(void)
{
   return(numberOfEvaluations);
}


// void setNumberOfVariables(int) method

/// This method sets a new number of variables in the objective function
///
/// @param newNumberOfVariables: New number of variables in the objective function.

void ObjectiveFunction::setNumberOfVariables(int newNumberOfVariables)
{
   if(newNumberOfVariables < 0)
   {
      std::cout << std::endl
                << "Error: ObjectiveFunction class. "
                << "void setNumberOfVariables(int) method."
                << std::endl
                << "Number of variables must be equal or greater than zero." 
                << std::endl
                << std::endl;

      exit(1);
   }

   numberOfVariables = newNumberOfVariables;

   // Initialize lower bound of domain to -1

   Vector<double> newLowerBound(numberOfVariables, -1.0);
   
   lowerBound = newLowerBound;

   // Initialize lower bound of domain to +1
   
   Vector<double> newUpperBound(numberOfVariables, 1.0);

   upperBound =  newUpperBound;
}


// void setLowerBound(Vector<double>)

/// This method sets a new lower bound in the objective function domain.
///
/// @param newLowerBound: New lower bound in the objective function domain.
///
/// @see setUpperBound(Vector<double>)
/// @see setDomain(Matrix<double>)

void ObjectiveFunction::setLowerBound(Vector<double> newLowerBound)
{
   if(newLowerBound.getSize() != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: ObjectiveFunction class. "
                << "void setLowerBound(Vector<double>) method."
                << std::endl
                << "Size must be equal to number of variables." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      // Check that lower bound of all variables is not greater than upper bound

      for(int i = 0; i < numberOfVariables; i++)
      {
         if(newLowerBound[i] > upperBound[i])
         {
            std::cout << std::endl
                      << "Error: ObjectiveFunction class. "
                      << "void setLowerBound(Vector<double>) method."
                      << std::endl
                      << "Lower bound of variable "<< i 
                      << " is greater than upper bound of that variable."
                      << std::endl << std::endl;

            exit(1);
         }
      }
   }
   
   // Set lower bound 

   lowerBound = newLowerBound;
}


// void setUpperBound(Vector<double>)

/// This method sets a new upper bound in the objective function domain.
///
/// @param newUpperBound: New upper bound in the objective function domain.
///
/// @see setLowerBound(Vector<double>)
/// @see setDomain(Matrix<double>)

void ObjectiveFunction::setUpperBound(Vector<double> newUpperBound)
{
   if(newUpperBound.getSize() != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: ObjectiveFunction class. "
                << "void setUpperBound(Vector<double>) method."
                << std::endl
                << "Size must be equal to number of variables." << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      // Check that upper bound of all variables is not less than lower bound

      for(int i = 0; i < numberOfVariables; i++)
      {
         if(newUpperBound[i] < lowerBound[i])
         {
            std::cout << std::endl
                      << "Error: ObjectiveFunction class. "
                      << "void setUpperBound(Vector<double>) method."
                      << std::endl
                      << "Upper bound of variable "<< i 
                      << " is less than lower bound of that variable."
                      << std::endl << std::endl;

            exit(1);
         }
      }
   }
   
   // Set upper bound 

   upperBound = newUpperBound;
}

// void setDomain(Matrix<double>)

/// This method sets a new objective function domain from a single Matrix.
/// Row 0 must contain the lower bound of the domain.
/// Row 1 must contain the upper bound of the domain.
///
/// @param newDomain: New objective function domain.
///
/// @see setLowerBound(Vector<double>)
/// @see setUpperBound(Vector<double>)

void ObjectiveFunction::setDomain(Matrix<double> newDomain)
{
   if(newDomain.getNumberOfRows() != 2)
   {
      std::cout << std::endl
                << "Error: ObjectiveFunction class. "
                << "void setDomain(Matrix<double>) method."
                << std::endl
                << "Number of rows must be 2." 
                << std::endl << std::endl;

      exit(1);
   }
   else if(newDomain.getNumberOfColumns() != numberOfVariables)
   {
      std::cout << std::endl
                << "Error: ObjectiveFunction class. "
                << "void setDomain(Matrix<double>) method."
                << std::endl
                << "Number of columns must be equal to number of variables." 
                << std::endl
                << std::endl;

      exit(1);
   }
   else
   {
      // Check that minimum of input variables is not greater than their maximum

      for(int i = 0; i < numberOfVariables; i++)
      {
         if(newDomain[0][i] > newDomain[1][i])
         {
            std::cout << std::endl
                      << "Error: ObjectiveFunction class. "
                      << "void setDomain(Matrix<double>) method."
                      << std::endl
                      << "Lower bound of input variable "<< i 
                      << " is greater than upper bound of that variable."
                      << std::endl << std::endl;

            exit(1);
         }
      }
   }
   
   // Set domain 

   for(int i = 0; i < numberOfVariables; i++)
   {
      lowerBound[i] = newDomain[0][i];
      upperBound[i] = newDomain[1][i];
   }   
}



// void setEpsilon(double) method

/// This method sets a new epsilon value to be used for numerical differentiation.
///
/// @param newEpsilon: New value for epsilon.
///
/// @see getGradient(Vector<double>)
/// @see getHessian(Vector<double>)

void ObjectiveFunction::setEpsilon(double newEpsilon)
{
   if(newEpsilon <= 0.0)
   {
      std::stringstream buffer;

      buffer << std::endl
             << "Error: ObjectiveFunction class. "
             << "void setEpsilon(double) method." << std::endl
             << "Epsilon must be greater than 0." << std::endl
             << std::endl;

      throw std::invalid_argument(buffer.str());
      //exit(1);
   }
   else
   {
      epsilon = newEpsilon;
   }
}


// void setNumberOfEvaluations(int) method

/// This method sets the number of calls to the getEvaluation(Vector<double>)
/// method to a new value. 
///
/// @param newNumberOfEvaluations: New number of calls
/// to the getEvaluation(Vector<double>) method.
///
/// @see getEvaluation(Vector<double>).

void ObjectiveFunction::setNumberOfEvaluations(int newNumberOfEvaluations)
{
   numberOfEvaluations = newNumberOfEvaluations;
}


// Vector<double> getGradient(void) method

/// This method returns the objective function gradient vector at a given argument.
/// It uses numerical differentiation.
///
/// @param argument: Point at which the gradient is to be computed.
///
/// @see getEvaluation(Vector<double>).
/// @see getHessian(Vector<double>).

Vector<double> ObjectiveFunction::getGradient(Vector<double> argument)
{
   Vector<double> gradient(numberOfVariables, 0.0);

   double evaluation1 = 0.0, evaluation2 = 0.0;

   for (int i = 0; i < numberOfVariables; i++)
   {
      // Add epsilon to argument

      argument[i] += epsilon;

      // Get evaluation

      evaluation2 = getEvaluation(argument);

      // Restart original argument

      argument[i] -= epsilon;

      // Substract epsilon from argument

      argument[i] -= epsilon;

      // Get evaluation

      evaluation1 = getEvaluation(argument);

      // Restart original argument

      argument[i] += epsilon;

      // Calculate derivative

      gradient[i] = (evaluation2 - evaluation1)/(2.0*epsilon);
   }

   return(gradient);
}


// double getGradientNorm(Vector<double>) method

/// This method returns the norm of the objective function gradient.
///
/// @param gradient: Objective function gradient Vector.
///
/// @see getGradient(Vector<double>).

double ObjectiveFunction::getGradientNorm(Vector<double> gradient)
{
   double gradientNorm = 0.0;

   int numberOfDimensions = gradient.getSize();

   double sum = 0.0;

   for(int i = 0; i < numberOfDimensions; i++)
   {
      sum += pow(gradient[i], 2);
   }

   gradientNorm = sqrt(sum);

   return(gradientNorm);
}


// Matrix<double> getHessian(Vector<double>)

/// This method returns the objective function Hessian matrix at a given argument.
/// It uses numerical differentiation.
///
/// @param argument: Point at which the Hessian is to be computed.
///
/// @see getEvaluation(Vector<double>).
/// @see getGradient(Vector<double>).

Matrix<double> ObjectiveFunction::getHessian(Vector<double> argument)
{
   Matrix<double> hessian(numberOfVariables, numberOfVariables, 0.0);

   double evaluation11 = 0.0, evaluation12 = 0.0;
   double evaluation21 = 0.0, evaluation22 = 0.0;

   // Obtain the upper part of the Hessian matrix

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = i; j < numberOfVariables; j++)
      {
         // Perturb argument components i and j

         argument[i] += epsilon;
         argument[j] += epsilon;

         // Calculate evaluation

         evaluation22 = getEvaluation(argument);

         // Restart argument components i and j

         argument[i] -= epsilon;
         argument[j] -= epsilon;


         // Perturb argument components i and j

         argument[i] += epsilon;
         argument[j] -= epsilon;

         // Calculate evaluation

         evaluation21 = getEvaluation(argument);

         // Restart argument components i and j

         argument[i] -= epsilon;
         argument[j] += epsilon;


         // Perturb argument components i and j

         argument[i] -= epsilon;
         argument[j] += epsilon;

         // Calculate evaluation

         evaluation12 = getEvaluation(argument);

         // Restart potential free parameters i and j

         argument[i] += epsilon;
         argument[j] -= epsilon;


         // Perturb potential free parameters i and j

         argument[i] -= epsilon;
         argument[j] -= epsilon;

         // Calculate evaluation

         evaluation11 = getEvaluation(argument);

         // Restart potential free parameters i and j

         argument[i] += epsilon;
         argument[j] += epsilon;

         // Calculate second derivative

         hessian[i][j]
         = (evaluation22 - evaluation21 - evaluation12 + evaluation11)/(4.0*pow(epsilon,2));
      }
   }

   // Obtain the rest of elements by symmetry

   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < i; j++)
      {
         hessian[i][j] = hessian[j][i];
      }
   }

   return(hessian);
}


// Matrix<double> getInverseHessian(Vector<double>) method

/// This method computes the Hessian at a given argument and then returns its
/// inverse.
///
/// @param argument: Point at which the inverse Hessian is to be computed.
///
/// @see getHessian(Vector<double>).

Matrix<double> ObjectiveFunction::getInverseHessian(Vector<double> argument)
{
   Matrix<double> inverseHessian(numberOfVariables, numberOfVariables, 0.0);
   
   Matrix<double> hessian = getHessian(argument);

   
   double hessianDeterminant = getDeterminant(hessian);

   if(hessianDeterminant == 0.0)
   {
      std::cout << "Error: ObjectiveFunction class. "
                << "Matrix<double> getInverseHessian(Vector<double>) method." 
                << std::endl
                << "Hessian matrix is singular." << std::endl
                << std::endl;
      
      exit(1);
   }
   
   // Get cofactor matrix
   
   Matrix<double> cofactor(numberOfVariables, numberOfVariables, 0.0);
                  
   Matrix<double> c(numberOfVariables-1, numberOfVariables-1, 0.0);

   for(int j = 0; j < numberOfVariables; j++) 
   {
      for (int i = 0; i < numberOfVariables; i++) 
      {
//         Form the adjoint a_ij
         int i1 = 0;

         for(int ii = 0; ii < numberOfVariables; ii++) 
         {
            if(ii == i)
            {
               continue;
            }
            
            int j1 = 0;

            for(int jj = 0; jj < numberOfVariables; jj++) 
            {
               if (jj == j)
               {
                  continue;
               }

               c[i1][j1] = hessian[ii][jj];
               j1++;
            }
            i1++;
         }

         double determinant = getDeterminant(c);

         cofactor[i][j] = pow(-1.0, i+j+2.0)*determinant;
      }
   }

   // Adjoint matrix is the transpose of cofactor matrix

   Matrix<double> adjoint(numberOfVariables, numberOfVariables, 0.0);
   
   double temp = 0.0;

   for(int i = 0; i < numberOfVariables; i++) 
   {
      for (int j = 0; j < numberOfVariables; j++) 
      {
         adjoint[i][j] = cofactor[j][i];
      }
   }

   // Inverse matrix is adjoint matrix divided by matrix determinant
   
   for(int i = 0; i < numberOfVariables; i++)
   {
      for(int j = 0; j < numberOfVariables; j++)
      {
         inverseHessian[i][j] = adjoint[i][j]/hessianDeterminant;
      }        
   } 
   
   
   return(inverseHessian);               
}


// void print(void) method

void ObjectiveFunction::print(Vector<double> argument)
{
   // Do nothing
}


// double getDeterminant(Matrix<double>) method

/// This method returns the determinant of a Matrix.
///
/// @param matrix: Matrix.

double ObjectiveFunction::getDeterminant(Matrix<double> matrix)
{
   double determinant = 0.0;

   int numberOfRows = matrix.getNumberOfRows();
   int numberOfColumns = matrix.getNumberOfColumns();

   if(numberOfRows != numberOfColumns)
   {
      std::cout << "Error: NewtonMethod class. "
                << "getDeterminant(Matrix<double>) method." << std::endl
                << "Matrix must be square" << std::endl
                << std::endl;
      
      exit(1);
   }

   if(numberOfRows == 0)
   {
      std::cout << "Error: NewtonMethod class. "
                << "getDeterminant(Matrix<double>) method." << std::endl
                << "Size of matrix is zero." << std::endl
                << std::endl;
      
      exit(1);                   
   }
   else if(numberOfRows == 1)
   {
//      std::cout << "Warning: NewtonMethod class. "
//                << "getDeterminant(Matrix<double>) method." << std::endl
//                << "Size of matrix is one." << std::endl;

      determinant = matrix[0][0];                   
   }
   else if(numberOfRows == 2)
   {
      determinant = matrix[0][0]*matrix[1][1] - matrix[1][0]*matrix[0][1];
   }
   else
   {
      for(int j1 = 0; j1 < numberOfRows; j1++) 
      {
         Matrix<double> subMatrix(numberOfRows-1, numberOfColumns-1, 0.0);     
     
         for(int i = 1; i < numberOfRows; i++) 
         {
            int j2 = 0;
      
            for (int j = 0; j < numberOfColumns; j++) 
            {
               if (j == j1)
               {
                  continue;
               }

               subMatrix[i-1][j2] = matrix[i][j];

               j2++;
            }
         }
   
         determinant += pow(-1.0, j1+2.0)*matrix[0][j1]*getDeterminant(subMatrix);    
      }
   }
      

   return(determinant);
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
