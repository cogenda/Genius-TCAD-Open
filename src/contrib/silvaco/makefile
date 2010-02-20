All: GENIUS

include ../../../make.defs

INCLUDES = -I.

GENIUS:   silvaco.o

silvaco.o       : silvaco.cc      silvaco.h

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


.PHONY : clean distclean

clean:
	@rm -f *.o


distclean:
	@rm -f *.o

