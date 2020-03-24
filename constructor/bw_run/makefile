CC = /opt/intel/bin/icpc
CC_FLAGS = -I/u/sciteam/skhan/eigen/ -I/opt/cray/hdf5/default/intel/15.0/include/ -L/opt/cray/hdf5/default/intel/15.0/lib/ -L/opt/intel/compilers_and_libraries/linux/lib/intel64_lin/ -lhdf5 -qopenmp -parallel  -mkl=parallel -std=c++11 -O3 -xhost -DUse_Indef -DUse_MKL -DNDEBUG -W

# File names
EXEC = run
SOURCES = $(wildcard ./src/*.cpp)
OBJECTS = $(SOURCES:.cpp=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(CC_FLAGS) -o $(EXEC)

# To obtain object files
%.o: %.cpp
	$(CC) -c $(CC_FLAGS) $< -o $@

# To remove generated files
clean:
	rm -f $(EXEC) $(OBJECTS)

