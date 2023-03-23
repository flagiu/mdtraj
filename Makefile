CXX = g++ -O3
SRC = src/mdtraj.cpp
OBJ = $(SRC:.cpp = .o)

mdtraj: $(OBJ)
	$(CXX) -o bin/mdtraj $(OBJ)

clean:
	rm -f core src/*.o bin/mdtraj
