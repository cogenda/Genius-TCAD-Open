include ../../../make.defs

INCLUDE = -I.

SRC     := csa.cc svd.c
OBJS   	:= csa.o  svd.o


.PHONY: clean clobber distclean

target := libcsa$(static_libext)


libcsa$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

svd.o : svd.c
	@echo "Compiling C (in "$(mode)" mode) "$<"..."
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

csa.o : csa.cc
	@echo "Compiling CXX (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

clean:
	@rm -f $(OBJS) libcsa.*

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libcsa.*

distclean:
	@$(MAKE) clobber

