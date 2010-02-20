include ../../../make.defs

DFISE_INCLUDE = -I.

SRCS    := $(wildcard *.cc)
OBJS    := $(patsubst %.cc, %$(obj-suffix), $(SRCS))


target := libdfise$(static_libext)


libdfise$(static_libext): $(OBJS)
	@$(AR) r $@ $(OBJS)

dfise.cc : dfise_lex.yy.c dfise_parser.tab.c
	@touch dfise.cc


dfise_lex.yy.c: dfise_lex.l
	@$(LEX) -l -i dfise_lex.l
	-@$(GENIUS_DIR)/bin/patch_lex lex.yy.c
	@mv lex.yy.c dfise_lex.yy.c

dfise_parser.tab.c : dfise_parser.y
	@$(YACC) -d dfise_parser.y
	@mv y.tab.h dfise_parser.tab.h
	@mv y.tab.c dfise_parser.tab.c

dfise_parser.tab.h : dfise_parser.tab.c

%$(obj-suffix) : %.cc
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(DFISE_INCLUD) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@

# generate c/c++ source depedence
.depends: $(SRCS)
	@if [ $(SYSTEM) != CYGWIN ]; then \
          $(CXX) -MM $(SRCS) $(DFISE_INCLUDE) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) > .depends ; \
        else \
          echo > .depends ; \
        fi

-include .depends

.PHONY: clean clobber distclean

clean:
	@rm -f $(OBJS) libdfise$(static_libext)

clobber:
	@$(MAKE) clean
	@rm -f *.o *.g.o *.pg.o *.sy.o
	@rm -f libdfise$(static_libext)
	@rm -f dfise_lex.yy.c  dfise_parser.tab.c dfise_parser.tab.h

distclean:
	@$(MAKE) clobber
