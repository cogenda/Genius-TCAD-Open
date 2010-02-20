include ../../../make.defs

SRCS    := qshep2d.f90 qshep3d.f90
OBJS    := qshep2d.o   qshep3d.o


target := libqshep$(static_libext)


libqshep$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

qshep2d.o: qshep2d.f90
	@echo "Compiling fortran (in "$(mode)" mode) "$<"..."
	@$(FC) -c -O2 $(DLLFLAG) $< -o $@

qshep3d.o: qshep3d.f90
	@echo "Compiling fortran (in "$(mode)" mode) "$<"..."
	@$(FC) -c -O2 $(DLLFLAG) $< -o $@

.PHONY: clean clobber distclean

clean:
	@rm -f $(OBJS) libqshep$(static_libext)

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libqshep$(static_libext)

distclean:
	@$(MAKE) clobber
