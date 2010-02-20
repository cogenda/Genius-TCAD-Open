ALL: target

include ../../../make.defs

SRCS    := $(wildcard *.cc)
OBJS    := $(patsubst %.cc, %$(obj-suffix), $(SRCS))


target : $(OBJS)


%$(obj-suffix) : %.cc
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE)  -c $< -o $@


# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends

.PHONY: clean clobber distclean

clean:
	@rm -f $(OBJS)

clobber:
	@$(MAKE) clean
	-@rm -f *.o *.g.o *.pg.o *.sy.o
	-@rm -f .depends

distclean:
	@$(MAKE) clobber
