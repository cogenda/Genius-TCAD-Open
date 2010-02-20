ALL: Genius

include ../../make.defs

SRCS   := physical_unit.cc  material_define.cc  material.cc
OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))
SUBDIRS := ${shell ls -F | grep /}

Genius: .depends adolc_init.o $(OBJS)
	@for dir in $(SUBDIRS) ; do \
		(cd $$dir ; make) || exit 1 ; \
	done

adolc_init.o : adolc_init.cc ../../include/math/adolc.h
	@echo "Compiling C++ (in "$(mode)" mode) with dynamic load support "$<"..."
	@$(CXX) $(CXXFLAGS) $(DLLFLAG) -c $< $(GENIUS_INCLUDE) $(PETSC_INCLUDE)

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@

# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends


clean:
	@cd Al;      make clean;
	@cd Ag;      make clean;
	@cd Au;      make clean;
	@cd Cu;      make clean;
	@cd ITO;     make clean;
	@cd TiSi2;   make clean;
	@cd PolySi;  make clean;
	@cd NPolySi; make clean;
	@cd PPolySi; make clean;
	@cd Elec;    make clean;
	@cd Si;      make clean;
	@cd GaAs;    make clean;
	@cd Ge;      make clean;
	@cd SiGe;    make clean;
	@cd AlGaAs;  make clean;
	@cd InGaAs;  make clean;
	@cd InP;     make clean;
	@cd InAs;    make clean;
	@cd InSb;    make clean;
	@cd InN;     make clean;
	@cd HgCdTe;  make clean;
	@cd 3C-SiC;  make clean;
	@cd 4H-SiC;  make clean;
	@cd S-SiO2;  make clean;
	@cd Air;     make clean;
	@cd SiO2;    make clean;
	@cd HfO2;    make clean;
	@cd Nitride; make clean;
	@cd Vacuum;  make clean;
	@cd PML;     make clean;
	@rm -f *.o
	@rm -f *.so

distclean:
	-@rm -f .depends
	@cd Al;      make distclean;
	@cd Ag;      make distclean;
	@cd Au;      make distclean;
	@cd Cu;      make distclean;
	@cd ITO;     make distclean;
	@cd TiSi2;   make distclean;
	@cd PolySi;  make distclean;
	@cd NPolySi; make distclean;
	@cd PPolySi; make distclean;
	@cd Elec;    make distclean;
	@cd Si;      make distclean;
	@cd GaAs;    make distclean;
	@cd Ge;      make distclean;
	@cd SiGe;    make distclean;
	@cd AlGaAs;  make distclean;
	@cd InGaAs;  make distclean;
	@cd InP;     make distclean;
	@cd InAs;    make distclean;
	@cd InSb;    make distclean;
	@cd InN;     make distclean;
	@cd HgCdTe;  make distclean;
	@cd 3C-SiC;  make distclean;
	@cd 4H-SiC;  make distclean;
	@cd S-SiO2;  make distclean;
	@cd Air;     make distclean;
	@cd SiO2;    make distclean;
	@cd HfO2;    make distclean;
	@cd Nitride; make distclean;
	@cd Vacuum;  make distclean;
	@cd PML;     make distclean;
