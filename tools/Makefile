CXX = g++
CXXFLAGS = -g -std=c++11 -Wall -Og -fsanitize=address
LIBS = -fsanitize=address
EIGENFLAGS = -I/usr/include/eigen3

all: permmatch

clean:
	rm permmatch permmatch.o pdb.o

permmatch: permmatch.o pdb.o
	$(CXX) $^ -o $@ $(LIBS)

%.o: %.cpp
	$(CXX) -c -o $@ $< $(CXXFLAGS) $(EIGENFLAGS)
