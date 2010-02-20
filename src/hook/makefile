ALL: Genius

include ../../make.defs

ifeq ($(SYSTEM),CYGWIN)
  HOOKS := shell_hook.o  rawfile_hook.o  gnuplot_hook.o  cv_hook.o  probe_hook.o  vtk_hook.o
  SRCS  := shell_hook.cc rawfile_hook.cc gnuplot_hook.cc cv_hook.cc probe_hook.cc vtk_hook.cc
else
  HOOKS := dlhook.o  shell_hook.so rawfile_hook.so gnuplot_hook.so cv_hook.so probe_hook.so vtk_hook.so
  SRCS  := dlhook.cc shell_hook.cc rawfile_hook.cc gnuplot_hook.cc cv_hook.cc probe_hook.cc vtk_hook.cc
endif


# we have to static link visit hook by now!
ifeq ($(enable-visit),yes)
    HOOKS := $(HOOKS)  visit_hook.o
    SRCS  := $(SRCS)   visit_hook.cc
endif

Genius: $(HOOKS)

# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) with dll support "$<"..."
	@$(CXX) $(CXXFLAGS) $(DLLFLAG) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@


shell_hook.so    :   shell_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

rawfile_hook.so  :   rawfile_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

gnuplot_hook.so  :   gnuplot_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

cv_hook.so       :   cv_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

probe_hook.so       :   probe_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

vtk_hook.so      :   vtk_hook.o
	@echo "Linking dynamic load library "$@"..."
	@$(CXX) $(CXXSHAREDFLAG) -o $@  $(DLLFLAG) $< $(LIBS)  $(CLIBS)

.PHONY : clean distclean

clean:
	-@rm -f *.o
	-@rm -f *.so

distclean:
	-@rm -f .depends
