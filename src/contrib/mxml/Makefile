include ../../../make.defs
ifeq ($(SYSTEM),CYGWIN)
    CFLAGS        += -DWIN32
    CXXFLAGS      += -DWIN32
endif

CFLAGS += -I.

SRC     := $(wildcard *.c)
OBJS	= mxml-attr.o mxml-entity.o mxml-file.o mxml-index.o \
	  mxml-node.o mxml-search.o mxml-set.o \
	  $(PUBLIBOBJS) mxml-private.o mxml-string.o

.PHONY: clean clobber distclean

%$(obj-suffix) : %.c config.h
	@echo "Compiling C (in "$(mode)" mode) "$<"..."
	@$(CC) $(CFLAGS) -c $< -o $@

%$(obj-suffix) : %.cpp config.h
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) -c $< -o $@

all : libmxml$(static_libext) MXMLUtil.o

libmxml$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

clean:
	@rm -f $(OBJS) libmxml.*

clobber:
	@$(MAKE) clean

distclean:
	@$(MAKE) clobber


