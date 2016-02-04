all: libPMIBenchmark.a

include ../make.defs

GENIUS_INCLUDE:=-I.. -I../include -I.

SRCS   := ../adolc_init.cc ../parser_parameter.cc ../PMI.cc  \
          PMI_benchmark_basic.cc \
          PMI_benchmark_band.cc \
          PMI_benchmark_mob.cc  \
          PMI_benchmark_thermal.cc \
          PMI_benchmark_avalanche.cc \
          PMI_benchmark_optical.cc

OBJS     := $(patsubst %.cc, %$(obj-suffix), $(SRCS))


.cc.o:
	@echo "Compiling C++ with dynamic load support "$<"..."
	@$(CXX) $(CXXFLAGS) $(GENIUS_INCLUDE) $(DLLFLAG) -c $< -o $@

#######################################################
libPMIBenchmark.so: $(OBJS)
	@echo "Linking PMI Benchmark files "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS)

libPMIBenchmark.a: $(OBJS)
	@echo "Linking PMI Benchmark files "$@
	@$(AR) rv $@ $(OBJS)

.PHONY: clean

clean:
	-@rm -f *.o *.so *.a

