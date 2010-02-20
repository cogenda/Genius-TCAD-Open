include ../../../make.defs

SRC     := $(wildcard *.cc)
OBJS   	:= $(patsubst %.cc, %$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

%$(obj-suffix) : %.cc
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) -c -I. $< -o $@

target := libexpreval$(static_libext)


libexpreval$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

clean:
	@rm -f $(OBJS) libexpreval.*
	@rm -f *~ *.o

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libexpreval.*

distclean:
	@$(MAKE) clobber

