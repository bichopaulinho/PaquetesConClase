/************************************************************************/
/*                                                                      */
/*   V E C T O R   C O N T A I N E R                                    */
/*                                                                      */
/*   Roberto Lopez                                                      */
/*   International Center for Numerical Methods in Engineering (CIMNE)  */
/*   Technical University of Catalonia (UPC)                            */
/*   Barcelona, Spain                                                   */
/*   E-mail: rlopez@cimne.upc.edu                                       */
/*                                                                      */
/************************************************************************/

#ifndef __VECTOR_H__
#define __VECTOR_H__

namespace Purple
{

/// This template class defines a vector for general purpose use.
///
/// @see Matrix.

template <class Type>
class Vector
{

private:

   /// Size of vector.

   int size;

   /// Pointer to a Type.

   Type* vector;

public:

   // CONSTRUCTORS

   Vector();

   Vector(int);

   Vector(int, const Type&);

   Vector(int, const Type*);

   Vector(const Vector&);


   // ASSINGMENT OPERATOR

   Vector& operator=(const Vector&);


   // REFERENCE OPERATORS

   inline Type& operator[](const int);

   inline const Type& operator[](const int) const;


   // METHODS

   inline int getSize();

   inline Type* begin();
   inline Type* end();


   // DESTRUCTOR

   ~Vector();
};


// CONSTRUCTORS


/// Default constructor. It creates a vector of size zero.

template <typename Type>
Vector<Type>::Vector() : size(0), vector(0)
{

}


/// Constructor. It creates a vector of size n, containing n copies of 
/// the default value for Type.
///
/// @param newSize: Size of Vector.

template <typename Type>
Vector<Type>::Vector(int newSize) : size(newSize), vector(new Type[newSize])
{

}


/// Constructor. It creates a vector of size n, containing n copies of 
/// the type value of Type. 
///
/// @param newSize: Size of Vector.
/// @param type: Value of Type.

template <typename Type>
Vector<Type>::Vector(int newSize, const Type& type) : size(newSize), vector(new Type[newSize])
{
   for(int i = 0; i < newSize; i++)
   {
      vector[i] = type;
   }
}


/// Constructor. It creates a vector of size n, containing n copies of 
/// the type value of Type. 
///
/// @param newSize: Size of Vector.
/// @param type: Value of Type.

template <typename Type>
Vector<Type>::Vector(int newSize, const Type* type) : size(newSize), vector(new Type[newSize])
{
   for(int i = 0; i < newSize; i++)
   {
      vector[i] = *type++;
   }
}


/// Copy constructor. It creates a copy of an existing Vector. 
///
/// @param oldVector: Vector to be copied.

template <typename Type>
Vector<Type>::Vector(const Vector<Type>& oldVector) : size(oldVector.size), vector(new Type[size])
{
   for(int i = 0; i < size; i++)
   {
      vector[i] = oldVector[i];
   }
}


// ASSIGNMENT OPERATORS

/// Assignment operator. It assigns to self a copy of an existing Vector.
///
/// @param oldVector: Vector to be assigned.

template <typename Type>
Vector<Type>& Vector<Type>::operator=(const Vector<Type>& oldVector)
{
   if (this != &oldVector)
   {
      if (size != oldVector.size)
      {
         if (vector != 0)
         {
            delete [] (vector);
         }

         size = oldVector.size;

         vector= new Type[size];
      }

      for (int i = 0; i < size; i++)
      {
         vector[i] = oldVector[i];
      }
   }

   return(*this);
}


// REFERENCE OPERATORS

/// Reference operator. 

template <typename Type>
inline Type& Vector<Type>::operator[](const int i) 
{
   return(vector[i]);
}

/// Reference operator. 

template <typename Type>
inline const Type& Vector<Type>::operator[](const int i) const 
{
   return(vector[i]);
}


// METHODS

// int getSize(void) method

/// This method returns the number of elements in the vector. 

template <typename Type>
inline int Vector<Type>::getSize()
{
   return(size);
}


// Type* begin(void) method

/// This method returns a pointer to the first element in the container.

template <typename Type>
inline Type* Vector<Type>::begin()
{
   return(vector);
}


// Type* end(void) method

/// This method returns a pointer to the last element in the container. 

template <typename Type>
inline Type* Vector<Type>::end()
{
   return(vector + size);
}


// DESTRUCTOR

/// Destructor. 

template <typename Type>
Vector<Type>::~Vector()
{
   if (vector != 0)
   {
      delete[] (vector);
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
