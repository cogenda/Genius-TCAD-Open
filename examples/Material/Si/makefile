all: libSi.so

include ../make.defs

GENIUS_INCLUDE:=-I../include -I../

SRCS   := $(wildcard *.cc)
SRCS   += ../adolc_init.cc ../PMI.cc

OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))

.cc.o:
	@echo "Compiling C++ with dynamic load support "$<"..."
	@$(CXX) $(CXXFLAGS) $(GENIUS_INCLUDE) $(DLLFLAG) -c $< -o $@  


libSi.so: $(OBJS)
	@echo "Linking material file "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS) 

.PHONY: clean

clean:
	-@rm -f *.o *.so

