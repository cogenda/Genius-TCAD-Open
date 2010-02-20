All: GENIUS

include ../../make.defs

GENIUS:   key.o pattern.o parser.o file_include.o pattern_xml.o

key.o          : key.cc       ../../include/parser/key.h
pattern.o      : pattern.cc   ../../include/parser/key.h ../../include/parser/pattern.h
parser.o       : parser.cc    ../../include/parser/key.h ../../include/parser/pattern.h ../../include/parser/parser.h
file_include.o : file_include.cc  ../../include/parser/file_include.h
pattern_xml.o : pattern_xml.cc ../../include/parser/pattern.h

.cc.o:
	@echo "Compiling C++ (in "$(mode)" mode) "$<"..."
	@$(CXX) $(CXXFLAGS) $(PETSC_INCLUDE) $(GENIUS_INCLUDE) -c $< -o $@

pattern.cc : pattern.yy.c pattern.tab.c
	@touch pattern.cc

parser.cc : input.yy.c input.tab.c
	@touch parser.cc

pattern.yy.c: pattern.l
	@$(LEX) -l -i pattern.l
	-@$(GENIUS_DIR)/bin/patch_lex lex.yy.c
	@mv lex.yy.c pattern.yy.c

pattern.tab.c : pattern.y
	@$(YACC) -d pattern.y
	@mv y.tab.h pattern.tab.h
	@mv y.tab.c pattern.tab.c

input.yy.c: input.l
	@$(LEX) -l -i input.l
	-@$(GENIUS_DIR)/bin/patch_lex lex.yy.c
	@mv lex.yy.c input.yy.c

input.tab.c : input.y
	@$(YACC) -d input.y
	@mv y.tab.h input.tab.h
	@mv y.tab.c input.tab.c





.PHONY : clean distclean

clean:
	@rm -f *.o


distclean:
#	@rm -f pattern.yy.c  pattern.tab.c pattern.tab.h input.yy.c input.tab.c input.tab.h

