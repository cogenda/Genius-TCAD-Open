include ../../make.defs


SRCS   := $(wildcard *.cc)
OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))

target := libgeom$(static_libext)

libgeom$(static_libext) : .depends $(OBJS)
	@$(AR) r $@ $(OBJS)

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
	-@rm -f libgeom$(static_libext) *.o

distclean:
	-@rm -f .depends






