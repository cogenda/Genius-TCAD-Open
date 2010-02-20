All:Genius

include ../make.defs
include ../objects

SRCS   := $(wildcard *.cc)

# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends


# note: For gcc, -Wl,--export-dynamic must be included in LDFLAGS for correctly loading hook plug-ins!
# Other compiler should also use equivalent flag!

Genius: .depends objs  genius$(EXE) libgenius$(static_libext)
#Genius: .depends objs  genius$(EXE)

libgenius$(shared_libext) : $(OBJ) $(LIB_OBJ)
	@echo "Create shared library "$@"..."
	-@$(LINK) $(CXXSHAREDFLAG) -o $@ $(OBJ)  $(LIBS) $(PETSC_LIB) $(CLIBS)

libgenius$(static_libext) : $(OBJ) $(LIB_OBJ)
	@echo "Create static library "$@"..."
	@$(AR) cr $@ $(OBJ)
	@../bin/concatlib $@ $(LIB_OBJ)

genius$(EXE):  main.o $(OBJ) $(LIB_OBJ)
	@echo "Linking executable binary "$@"..."
	@$(LINK) $(CXXFLAGS) $(LDFLAGS) -o genius$(EXE)  main.o $(OBJ)  $(LIBS) $(PETSC_LIB) $(CLIBS)
	@if test $(METHOD) = 'opt' ; then \
	  $(STRIP) genius$(EXE); \
	fi

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@

objs:
	@cd contrib;            make;
	@cd base;               make;
	@cd fe;                 make;
	@cd geom;               make;
	@cd interpolation;      make;
	@cd math;               make;
	@cd mesh;               make;
	@cd meshgen;            make;
	@cd material;           make;
	@cd parser;             make;
	@cd partition;          make;
	@cd quadrature;         make;
	@cd solution;           make;
	@cd solver;             make;
	@cd hook;               make;
	@cd source;             make;
	@cd utils;              make;


.PHONY : clean distclean

clean:
	@cd contrib;            make clean;
	@cd base;               make clean;
	@cd fe;                 make clean;
	@cd geom;               make clean;
	@cd interpolation;      make clean;
	@cd math;               make clean;
	@cd mesh;               make clean;
	@cd meshgen;            make clean;
	@cd material;           make clean;
	@cd parser;             make clean;
	@cd partition;          make clean;
	@cd quadrature;         make clean;
	@cd solution;           make clean;
	@cd solver;             make clean;
	@cd hook;               make clean;
	@cd source;             make clean;
	@cd utils;              make clean;
	@rm -f *.o
	@rm -f genius$(EXE)


distclean:
	@cd contrib;            make distclean;
	@cd base;               make distclean;
	@cd fe;                 make distclean;
	@cd geom;               make distclean;
	@cd interpolation;      make distclean;
	@cd math;               make distclean;
	@cd mesh;               make distclean;
	@cd meshgen;            make distclean;
	@cd material;           make distclean;
	@cd parser;             make distclean;
	@cd partition;          make distclean;
	@cd quadrature;         make distclean;
	@cd solution;           make distclean;
	@cd solver;             make distclean;
	@cd hook;               make distclean;
	@cd source;             make distclean;
	@cd utils;              make distclean;
	-@rm -f .depends


