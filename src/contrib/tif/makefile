All: GENIUS

include ../../../make.defs

INCLUDES = -I.

GENIUS:   tif_tool.o

tif_tool.o       : tif_tool.cc      tif_data.h

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(INCLUDES) $(GENIUS_INCLUDE) -c $< -o $@

tif_tool.cc : tif_lex.yy.c tif_parser.tab.c
	@touch tif_tool.cc


tif_lex.yy.c: tif_lex.l
	@$(LEX) -l -i tif_lex.l
	-@$(GENIUS_DIR)/bin/patch_lex lex.yy.c
	@mv lex.yy.c tif_lex.yy.c

tif_parser.tab.c : tif_parser.y
	@$(YACC) -d tif_parser.y
	@mv y.tab.h tif_parser.tab.h
	@mv y.tab.c tif_parser.tab.c


.PHONY : clean distclean

clean:
	@rm -f *.o


distclean:
	@rm -f tif_lex.yy.c  tif_parser.tab.c tif_parser.tab.h

