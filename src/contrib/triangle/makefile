include ../../../make.defs 

target := libctri$(static_libext)

libtriangle$(static_libext): triangle.o
	        @$(AR) r $@ triangle.o

triangle.o  : triangle.c
	@echo "Compiling C "$<"..."
	@$(CC)  -O0 -c -DTRILIBRARY $(DLLFLAG) triangle.c -I.


clean:
	-@rm -f *.o


distclean:
	-@rm -f *.o
