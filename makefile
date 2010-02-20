all:
	@cd src;    make;

include make.defs



install: src/genius$(EXE)
	@cp src/genius$(EXE) bin/genius$(EXE)
	@cp src/material/*.so lib/
	@if [ "$(SYSTEM)" != "CYGWIN" ]; then \
          cp src/hook/*.so lib/;  \
        fi;


.PHONY : clean distclean


clean:
	@cd src;      make clean;
	-@rm -f bin/genius.$(SYSTEM)$(EXE)
	-@rm -f lib/*.so


distclean:
	@cd src;      make distclean;
	-@rm -f config.h
	-@rm -f make.defs
	-@rm -f config.log
	-@rm -f config.status

