CXX = g++
# -O2 may trigger optimization bug (?)
CXXFLAGS = -g -std=c++11 -Wall -Og -fsanitize=address
LIBS = -fsanitize=address
EIGENFLAGS = -I/usr/include/eigen3

nucfepgen: nucfepgen.o pdb.o topology.o assign.o select.o
	$(CXX) -o $@ $^ $(LIBS)

nucfepgen.o: topology.hpp pdb.hpp assign.hpp select.hpp

clean:
	rm -f *.o nucfepgen

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(EIGENFLAGS) -c -o $@ $<

