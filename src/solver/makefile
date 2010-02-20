ALL: Genius

include ../../make.defs

SRCS   := $(wildcard *.cc)
OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))

# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends


Genius: .depends $(OBJS)
	@cd doping_analytic; make;
	@cd mole_analytic;   make;
	@cd poisson;         make;
	@cd stress_solver;   make;
	@cd ddm_common;      make;
	@cd ddm1;            make;
	@cd ddm2;            make;
	@cd ebm3;            make;
	@cd hall;            make;
	@cd mix_common;      make;
	@cd mix1;            make;
	@cd mix2;            make;
	@cd mix3;            make;
	@cd mixA1;           make;
	@cd mixA2;           make;
	@cd mixA3;           make;
	@cd ddm_ac;          make;
	@cd emfem2d;         make;
	@cd ray_tracing;     make;


.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@

.PHONY : clean distclean

clean:
	-@rm -f *.o
	@cd doping_analytic; make clean;
	@cd mole_analytic;   make clean;
	@cd poisson;         make clean;
	@cd stress_solver;   make clean;
	@cd ddm_common;      make clean;
	@cd ddm1;            make clean;
	@cd ddm2;            make clean;
	@cd ebm3;            make clean;
	@cd hall;            make clean;
	@cd mix_common;      make clean;
	@cd mix1;            make clean;
	@cd mix2;            make clean;
	@cd mix3;            make clean;
	@cd mixA1;           make clean;
	@cd mixA2;           make clean;
	@cd mixA3;           make clean;
	@cd ddm_ac;          make clean;
	@cd emfem2d;         make clean;
	@cd ray_tracing;     make clean;

distclean:
	-@rm -f .depends
	@cd doping_analytic; make distclean;
	@cd mole_analytic;   make distclean;
	@cd poisson;         make distclean;
	@cd stress_solver;   make distclean;
	@cd ddm_common;      make distclean;
	@cd ddm1;            make distclean;
	@cd ddm2;            make distclean;
	@cd ebm3;            make distclean;
	@cd hall;            make distclean;
	@cd mix_common;      make distclean;
	@cd mix1;            make distclean;
	@cd mix2;            make distclean;
	@cd mix3;            make distclean;
	@cd mixA1;           make distclean;
	@cd mixA2;           make distclean;
	@cd mixA3;           make distclean;
	@cd ddm_ac;          make distclean;
	@cd emfem2d;         make distclean;
	@cd ray_tracing;     make distclean;
