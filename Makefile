### Fortran Compiler
FC = mpif90
#FC = ftn

### Compilation flags
FFLAGS := -O3 -ffixed-line-length-none -fdefault-real-8 

BIG  := #-mcmodel=medium
DBG  := #-g #-traceback
PROF := #-pg
OMP  := -fopenmp

### Generate script-dependent target (exe)
num := $(shell grep -o 'case_num = .*' param.f90 | sed -n 's/case_num = //p' | tr -d '"/')
TARGET := SLD$(num)

### Source scripts
SRC = param.f90 common.f90 zero.f90 core.f90 output.f90 main.f90

OBJ = $(SRC:.f90=.o)

all: $(TARGET)

$(TARGET): $(OBJ)
	$(FC) $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) -o $@ $(OBJ) $(LIB)

.PHONY: clean  
clean:
	rm -f *.o *.mod $(TARGET) *~ .depend

%.o: %.f90
	$(FC) $(INCLUDE) -c -o $@ $(FFLAGS) $(BIG) $(DBG) $(PROF) $(OMP) $<

.depend dep:
	./.makedepo $(SRC) > .depend

include .depend
