CC = /opt/intel/bin/icpc
CC_FLAGS = -I/media/bkcgroup/Install/eigen-eigen-dc6cfdf9bcec/ -I/usr/include/hdf5/serial/ -L/usr/lib/x86_64-linux-gnu/hdf5/serial/ -lhdf5 -qopenmp -parallel  -mkl=parallel -std=c++11 -O3 -xhost -DUse_Indef -DUse_MKL -DNDEBUG -W 

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
