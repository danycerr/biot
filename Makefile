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

clean:
	$(RM) $(OBJS)

resu_clean:
	$(RM) $(RESU_OBJS)
	
depend: .depend

.depend: $(SRCS)
	rm -f ./.depend
	$(CXX) $(CXXFLAGS) -MM $^ > ./.depend;

include .depend
	
.PHONY: all clean resu_clean depend
