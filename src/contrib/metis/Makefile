include ../../../make.defs

INCLUDE = -I. 

SRC     := $(wildcard *.c)
OBJS   	:= $(patsubst %.c, %$(obj-suffix), $(SRC))


.PHONY: clean clobber distclean

%$(obj-suffix) : %.c
	@echo "Compiling C (in "$(mode)" mode) "$<"..."
	@$(CC) $(CFLAGS) $(INCLUDE) -c $< -o $@

target := libmetis$(static_libext)


libmetis$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

clean:
	@rm -f $(OBJS) libmetis.*

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libmetis.*

distclean:
	@$(MAKE) clobber

