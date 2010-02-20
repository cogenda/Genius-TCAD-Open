ALL: brkpnts.h

include ../../../make.defs


brkpnts.h: brkpnts$(EXE)
	@echo "Execute "$<"..."
	@./brkpnts$(EXE)
	@cp brkpnts.h ../../../include/math

brkpnts$(EXE): brkpnts.c
	@echo "Compiling "$<"..."
	@$(CC) -o brkpnts$(EXE) brkpnts.c $(MATHLIB)

.PHONY : clean distclean

clean:
	@rm -f *.o brkpnts$(EXE) brkpnts.txt brkpnts.h

distclean:
	@rm -f *.o brkpnts$(EXE) brkpnts.txt brkpnts.h
	@rm -f ../../../include/math/brkpnts.h
