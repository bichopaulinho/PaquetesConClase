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
OUTPUTFILE  = ./lib/libRandomSearch.a

objects = ./Purple/bin/ObjectiveFunction.o ./Purple/bin/OptimizationAlgorithm.o ./Purple/bin/RandomSearch.o ./Purple/bin/OptimizationAlgorithm/GradientDescent.o ./Purple/bin/ObjectiveFunction/DeJongFunction.o
CXXFLAGS = -fPIC -Wall
CXX=g++
# Default target
.PHONY: all
all: $(OUTPUTFILE)

$(OUTPUTFILE): $(objects)
	ar rcs $@ $^

 
./Purple/bin/ObjectiveFunction.o: ./Purple/ObjectiveFunction/ObjectiveFunction.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/ObjectiveFunction/ObjectiveFunction.cpp -o  ./Purple/bin/ObjectiveFunction.o


./Purple/bin/OptimizationAlgorithm.o: ./Purple/OptimizationAlgorithm/OptimizationAlgorithm.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/OptimizationAlgorithm.cpp -o  ./Purple/bin/OptimizationAlgorithm.o

./Purple/bin/RandomSearch.o: ./Purple/OptimizationAlgorithm/RandomSearch.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/RandomSearch.cpp -o ./Purple/bin/RandomSearch.o

./Purple/bin/OptimizationAlgorithm/GradientDescent.o: ./Purple/OptimizationAlgorithm/GradientDescent.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/GradientDescent.cpp -o ./Purple/bin/OptimizationAlgorithm/GradientDescent.o

./Purple/bin/ObjectiveFunction/DeJongFunction.o: ./Purple/ObjectiveFunction/DeJongFunction.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/ObjectiveFunction/DeJongFunction.cpp -o ./Purple/bin/ObjectiveFunction/DeJongFunction.o

# files is required; this is handled by make's database of
# implicit rules

.PHONY: clean 
clean:
	for file in $(CLEANEXTS); do rm -f ./lib/*.$$file; rm -f ./Purple/bin/*.$$file; done





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
