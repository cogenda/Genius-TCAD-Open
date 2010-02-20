include ../../../make.defs

INCLUDES = -I. 

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

%$(obj-suffix) : %.c
	@echo "Compiling C (in "$(mode)" mode) "$<"..."
	@$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

target := libmetis$(static_libext)


libcuba$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

clean:
	@rm -f $(OBJS) libcuba.*

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libcuba.*

distclean:
	@$(MAKE) clobber

