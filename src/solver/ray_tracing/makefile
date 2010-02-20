ALL: Genius

include ../../../make.defs

SRCS   := $(wildcard *.cc)
OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))


Genius : .depends $(OBJS)

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


.PHONY : clean distclean

clean:
	-@rm -f *.o

distclean:
	-@rm -f .depends

