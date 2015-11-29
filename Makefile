################################################################################
#                                                                              # 
#   Purple: An Open Source Numerical Optimization C++ Library                  #
#                                                                              #
#   R A N D O M   S E A R C H   M A K E F I L E                                #
#                                                                              #
#   Roberto Lopez                                                              # 
#   International Center for Numerical Methods in Engineering (CIMNE)          #
#   Technical University of Catalonia (UPC)                                    #
#   Barcelona, Spain                                                           #
#   E-mail: rlopez@cimne.upc.edu                                               #
#                                                                              #
################################################################################



# Como hacer el makefile: https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch01s17.html
# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a 

# Specify the target file and the install directory
OUTPUTFILE  = libRandomSearch.a

objects = ObjectiveFunction.o OptimizationAlgorithm.o RandomSearch.o
CXXFLAGS = -fPIC -Wall
CXX=g++
# Default target
.PHONY: all
all: $(OUTPUTFILE)

$(OUTPUTFILE): $(objects)
	ar rcs $@ $^

 
ObjectiveFunction.o: ./Purple/ObjectiveFunction/ObjectiveFunction.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/ObjectiveFunction/ObjectiveFunction.cpp -o  ObjectiveFunction.o


OptimizationAlgorithm.o: ./Purple/OptimizationAlgorithm/OptimizationAlgorithm.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/OptimizationAlgorithm.cpp -o  OptimizationAlgorithm.o

RandomSearch.o: ./Purple/OptimizationAlgorithm/RandomSearch.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/RandomSearch.cpp -o RandomSearch.o

# files is required; this is handled by make's database of
# implicit rules

.PHONY: clean 
clean:
	for file in $(CLEANEXTS); do rm -f *.$$file; done





# Objective function



# Purple: An Open Source Numerical Optimization C++ Library.
# Copyright (C) 2006 Roberto Lopez
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
