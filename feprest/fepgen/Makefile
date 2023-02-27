CXX = g++
# Known issue: old GCC in Ubuntu 14 with -O2 or more may trigger optimization bug (?). Recent GCCs do not have the issue.
# I am using -fsanitize=address to catch in case there should exist out-of-bounds bug issue, but it's not mandatory.
# If you are using compilers other than g++ you can ignore these LIBS or CXXFLAGS.
CXXFLAGS = -g -std=c++11 -Wall -Og -fsanitize=address
LIBS = -fsanitize=address
EIGENFLAGS = -I/usr/include/eigen3

fepgen: fepgen.o pdb.o topology.o assign.o select.o
	$(CXX) -o $@ $^ $(LIBS)

fepgen.o: topology.hpp pdb.hpp assign.hpp select.hpp

clean:
	rm -f *.o fepgen

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(EIGENFLAGS) -c -o $@ $<

