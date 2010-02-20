ALL: ../libGe.so

include ../../../make.defs

SRCS   := $(wildcard *.cc)
SRCS   += ../adolc_init.cc
SRCS   += ../PMI.cc

OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))

../libGe.so:  libGe.so
	@cp $< ..

libGe.so: .depends $(OBJS)
	@echo "Linking material file "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS) $(MATHLIB)

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) with dynamic load support "$<"..."
	@$(CXX) $(CXXFLAGS) $(DLLFLAG) -c $< -o $@ $(GENIUS_INCLUDE) $(PETSC_INCLUDE)


# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends


.PHONY : clean distclean

clean:
	-@rm -f *.o *.so

distclean:
	-@rm -f .depends

