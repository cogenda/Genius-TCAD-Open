include ../../../make.defs

INCLUDE = -I. 

SRC     := $(wildcard *.cpp)
OBJS   	:= $(patsubst %.cpp, %$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

%$(obj-suffix) : %.cpp
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CFLAGS) $(INCLUDE) -c $< -o $@

target := libann$(static_libext)


libann$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

clean:
	@rm -f $(OBJS) libann.*

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libann.*

distclean:
	@$(MAKE) clobber

