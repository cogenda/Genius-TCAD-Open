all: libPMIBenchmarkMob.so  libPMIBenchmarkMob.a libPMIBenchmarkOptical.so libPMIBenchmarkOptical.a

include make.defs

GENIUS_INCLUDE:=-I./include -I.

SRCS   := adolc_init.cc parser_parameter.cc PMI.cc PMI_benchmark_mob.cc PMI_benchmark_optical.cc

OBJS   := $(patsubst %.cc, %$(obj-suffix), $(SRCS))

.cc.o:
	@echo "Compiling C++ with dynamic load support "$<"..."
	@$(CXX) $(CXXFLAGS) $(GENIUS_INCLUDE) $(DLLFLAG) -c $< -o $@


libPMIBenchmarkMob.so: $(OBJS)
	@echo "Linking PMI Benchmark Mobility file "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS)

libPMIBenchmarkMob.a: $(OBJS)
	@echo "Linking PMI Benchmark Mobility file "$@
	@$(AR) rv $@ $(OBJS)

libPMIBenchmarkOptical.so: $(OBJS)
	@echo "Linking PMI Benchmark Optical file "$@
	@$(CXX) $(CXXSHAREDFLAG) -o $@ $(OBJS)

libPMIBenchmarkOptical.a: $(OBJS)
	@echo "Linking PMI Benchmark Optical file "$@
	@$(AR) rv $@ $(OBJS)

.PHONY: clean

clean:
	-@rm -f *.o *.so *.a

