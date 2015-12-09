/************************************************************************/
/*                                                                      */
/*   M A T R I X   C O N T A I N E R                                    */
/*                                                                      */
/*   Roberto Lopez                                                      */
/*   International Center for Numerical Methods in Engineering (CIMNE)  */
/*   Technical University of Catalonia (UPC)                            */
/*   Barcelona, Spain                                                   */
/*   E-mail: rlopez@cimne.upc.edu                                       */
/*                                                                      */
/************************************************************************/


#ifndef __MATRIX_H__
#define __MATRIX_H__

namespace Purple
{

/// This template class defines a matrix for general purpose use.
///
/// @see Vector.

template <class Type>
class Matrix 
{

private:

   /// Number of rows in matrix.

   int numberOfRows;

   /// Number of columns in matrix.

   int numberOfColumns;

   /// Double pointer to a Type.

   Type** matrix;

public:

   // CONSTRUCTORS

   Matrix();

   Matrix(int, int);

   Matrix(int, int, const Type&);

   Matrix(int, int, const Type*);

   Matrix(const Matrix&);


   // ASSIGNMENT OPERATOR

   Matrix& operator=(const Matrix&);


   // REFERENCE OPERATORS

   inline Type* operator[](const int);

   inline const Type* operator[](const int) const;


   // METHODS

   inline int getNumberOfRows() const;
   inline int getNumberOfColumns() const;


   // DESTRUCTOR

   ~Matrix();
};


// CONSTRUCTORS


/// Default constructor. It creates a matrix with zero rows and zero 
/// columns.

template <class Type>
Matrix<Type>::Matrix() : numberOfRows(0), numberOfColumns(0), matrix(0) 
{

}


/// Constructor. It creates a matrix with n rows and m columns, containing 
/// n*m copies of the default value for Type. 
///
/// @param newNumberOfRows: Number of rows in Matrix.
/// @param newNumberOfColumns: Number of columns in Matrix.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), matrix(new Type*[newNumberOfRows])
{
   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for (int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }
}


/// Constructor. It creates a matrix with n rows and m columns, containing n*m copies of 
/// the type value of Type. 
///
/// @param newNumberOfRows: Number of rows in Matrix.
/// @param newNumberOfColumns: Number of columns in Matrix.
/// @param type: Value of Type.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns, const Type& type) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), matrix(new Type*[newNumberOfRows])
{
   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for (int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   for (int i = 0; i < numberOfRows; i++)
   {
      for (int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = type;
      }
   }
}


/// Constructor. It creates a matrix with n rows and m columns, containing n*m copies of 
/// the type value of Type. 
///
/// @param newNumberOfRows: Number of rows in Matrix.
/// @param newNumberOfColumns: Number of columns in Matrix.
/// @param type: Value of Type.

template <class Type>
Matrix<Type>::Matrix(int newNumberOfRows, int newNumberOfColumns, const Type* type) 
: numberOfRows(newNumberOfRows), numberOfColumns(newNumberOfColumns), matrix(new Type*[newNumberOfRows])
{
   // Construct matrix

   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for (int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   // Set all elements of matrix to the type value of Type

   for (int i = 0; i < numberOfRows; i++)
   {
      for (int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = *type++;
     }
   }
}


/// Copy constructor. It creates a copy of an existing Matrix. 
///
/// @param oldMatrix: Matrix to be copied.

template <class Type>
Matrix<Type>::Matrix(const Matrix& oldMatrix) 
: numberOfRows(oldMatrix.numberOfRows), numberOfColumns(oldMatrix.numberOfColumns), matrix(new Type*[numberOfRows])
{
   // Construct matrix

   matrix[0] = new Type[numberOfColumns*numberOfRows];

   for (int i = 1; i < numberOfRows; i++)
   {
      matrix[i] = matrix[i-1] + numberOfColumns;
   }

   // Set all elements of matrix to the old matrix type values

   for (int i = 0; i < numberOfRows; i++)
   {
      for (int j = 0; j < numberOfColumns; j++)
      {
         matrix[i][j] = oldMatrix[i][j];
      }
   }
}


// ASSIGNMENT OPERATORS

/// Assignment operator. It assigns to self a copy of an existing Matrix.
///
/// @param oldMatrix: Matrix to be assigned.

template <class Type>
Matrix<Type>& Matrix<Type>::operator=(const Matrix<Type>& oldMatrix)
{
   if (this != &oldMatrix) 
   {
      if (numberOfRows != oldMatrix.numberOfRows || numberOfColumns != oldMatrix.numberOfColumns) 
      {
         if (matrix != 0) 
         {
            delete[] (matrix[0]);

            delete[] (matrix);
         }

         numberOfRows = oldMatrix.numberOfRows;

         numberOfColumns = oldMatrix.numberOfColumns;

         matrix = new Type*[numberOfRows];

         matrix[0] = new Type[numberOfColumns*numberOfRows];
      }

      for (int i = 1; i < numberOfRows; i++)
      {
         matrix[i] = matrix[i-1] + numberOfColumns;
      }

      // Set all elements of matrix to the old matrix type values

      for (int i = 0; i < numberOfRows; i++)
      {
         for (int j = 0; j < numberOfColumns; j++)
         {
            matrix[i][j] = oldMatrix[i][j];
         }
      }
   }

   return(*this);
}


// REFERENCE OPERATORS

/// Reference operator.  

template <class Type>
inline Type* Matrix<Type>::operator[](const int i) 
{
   return(matrix[i]);
}


/// Reference operator.  

template <class Type>
inline const Type* Matrix<Type>::operator[](const int i) const
{
   return(matrix[i]);
}


// METHODS

// int getNumberOfRows(void) method

/// This method returns the number of rows in the matrix. 

template <class Type>
inline int Matrix<Type>::getNumberOfRows() const
{
   return(numberOfRows);
}


// int getNumberOfColumns(void) method

/// This method returns the number of columns in the matrix. 

template <class Type>
inline int Matrix<Type>::getNumberOfColumns() const
{
   return(numberOfColumns);
}


// DESTRUCTOR

/// Destructor. 

template <class Type>
Matrix<Type>::~Matrix()
{
   if (matrix != 0) 
   {
      delete[] (matrix[0]);

      delete[] (matrix);
   }
}

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
