# CXXFLAGS
INCLUDES =-I/u/archive/agip/cerroni/software/eigen/eigen/
CXX=g++
SRC_pol = polaretto.cpp
OBJ_pol = $(SRC_pol:.cpp=.o)

polaretto: $(OBJ_pol)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ_pol) 

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<