CC = icpc
CC_FLAGS = -I/home/aakhan3/eigen -I/usr/include/hdf5 -lhdf5 -qopenmp -parallel  -mkl=parallel -std=c++11 -O3 -DUse_Indef -DUse_MKL -DNDEBUG -W

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

