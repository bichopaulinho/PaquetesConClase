# Compilación de la versión ampliada de la biblioteca Purple en estático

# Como hacer el makefile: https://www.safaribooksonline.com/library/view/c-cookbook/0596007612/ch01s17.html
# Specify extensions of files to delete when cleaning
CLEANEXTS   = o a 

# Specify the target file and the install directory
OUTPUTFILE  = ./lib/libPurple.a

objects = ./Purple/bin/ObjectiveFunction.o ./Purple/bin/OptimizationAlgorithm.o ./Purple/bin/GradientDescent.o ./Purple/bin/SumaCuadrados.o
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

./Purple/bin/GradientDescent.o: ./Purple/OptimizationAlgorithm/GradientDescent.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/OptimizationAlgorithm/GradientDescent.cpp -o ./Purple/bin/GradientDescent.o

./Purple/bin/DeJongFunction.o: ./Purple/ObjectiveFunction/DeJongFunction.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/ObjectiveFunction/DeJongFunction.cpp -o ./Purple/bin/DeJongFunction.o

./Purple/bin/SumaCuadrados.o: ./Purple/ObjectiveFunction/SumaCuadrados.cpp
	$(CXX) $(CXXFLAGS) -c ./Purple/ObjectiveFunction/SumaCuadrados.cpp -o ./Purple/bin/SumaCuadrados.o

# files is required; this is handled by make's database of
# implicit rules

.PHONY: clean 
clean:
	for file in $(CLEANEXTS); do rm -f ./lib/*.$$file; rm -f ./Purple/bin/*.$$file; done

