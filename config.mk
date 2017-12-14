# getfem
CXXFLAGS+=$(shell getfem-config --cflags)
LDFLAGS+=$(shell getfem-config --libs)

# superlu
CXXFLAGS+=-DGMM_USES_SUPERLU -I$(mkSuperluInc)
LDFLAGS+=-L$(mkSuperluLib) -lsuperlu
LDFLAGS+=-L$(mkQhullLib)

CXXFLAGS+=-std=c++14
