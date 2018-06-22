# getfem
CXXFLAGS+=$(shell getfem-config --cflags)
LDFLAGS+=$(shell getfem-config --libs)

# superlu
CXXFLAGS+=-DGMM_USES_SUPERLU -I$(mkSuperluInc)
LDFLAGS+=-L$(mkSuperluLib) -lsuperlu
LDFLAGS+=-L$(mkQhullLib)

CXXFLAGS+=-std=c++14


LDFLAGS += -L/opt/lib/samg/ -lamg -liomp5

CXXFLAGS += -I ${SAMG}/
CXXFLAGS+= -DSAMG_UNIX_LINUX -DSAMG_LCASE_USCORE -DPYRAMID_TRIANGULAR_FACETS
