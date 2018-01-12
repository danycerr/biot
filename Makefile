include config.mk

SRCS=$(wildcard *.cpp)
HDRS=$(wildcard *.hpp)
OBJS=$(SRCS:.cpp=.o)
FIN_EXEC=lev_set

all: $(FIN_EXEC)

$(FIN_EXEC): $(OBJS)
	$(CXX) -o $@ $(OBJS) $(LDFLAGS) $(LIBS)
        
%.o: %.cpp config.mk
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c -o $@ $<

%.cpp: %.hpp
	touch $@
	
main.cpp: biot_ls.hpp
	touch $@

clean:
	$(RM) $(OBJS)

resu_clean:
	$(RM) $(RESU_OBJS)
	
.PHONY: all clean resu_clean 
