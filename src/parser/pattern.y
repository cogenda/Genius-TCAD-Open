%{

#include <stdio.h>
#include <string.h>

extern int yylineno;
extern int yylex();
static int yyerror(void * dummy, const char *s);

static struct Parser::PatternCard card;
static Parser::Parameter   parameter;

//#define VERBOSE

%}

%start inputfile

%union{
    bool   bval;
    int    ival;
    double dval;
    char   sval[80];
}


%parse-param { void * dummy }

%token <ival> INTEGER_VALUE
%token <dval> REAL_VALUE
%token <sval> STRING_VALUE
%token <bval> BOOL_VALUE
%token <ival> INTEGER REAL STRING BOOL ENUM
%token <ival> BLANK COMMENT BAD_WORD



%%

inputfile
    :   command
    |   inputfile command
    ;

command :  STRING_VALUE '{' parameters '}'
{
#ifdef VERBOSE
	printf("line %d: -- YACC command %s --\n",yylineno-1,$1);
#endif

	//change key to upper case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(islower($1[c]))
            $1[c]=toupper($1[c]);
        }

	card._key = $1;
	Parser::Pattern * pattern = (Parser::Pattern *)dummy;
	pattern->_pattern_card_map.insert(std::pair<const std::string, const Parser::PatternCard>($1,card));
	card.clear();
}

	|  STRING_VALUE '{'  '}'
{
#ifdef VERBOSE
	printf("line %d: -- YACC command %s --\n",yylineno-1,$1);
#endif

	//change key to upper case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(islower($1[c]))
            $1[c]=toupper($1[c]);
        }

	card._key = $1;
	Parser::Pattern * pattern = (Parser::Pattern *)dummy;
	pattern->_pattern_card_map.insert(std::pair<const std::string, const Parser::PatternCard>($1,card));
	card.clear();
}
	|  COMMENT
        ;


parameters : parameter
	   | parameters parameter
	   ;

parameter  : STRING_VALUE   BOOL     BOOL_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC patameter %s BOOL %d --\n",yylineno-1,$1,int($3));
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	parameter.set_name($1);
	parameter.set_bool($3);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}
	   | STRING_VALUE   INTEGER  INTEGER_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC patameter %s INT %d --\n",yylineno-1,$1,int($3));
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	parameter.set_name($1);
	parameter.set_int($3);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}
	   | STRING_VALUE   REAL     REAL_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC patameter %s REAL %e --\n",yylineno-1,$1,$3);
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	parameter.set_name($1);
	parameter.set_real($3);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}
           | STRING_VALUE    REAL     INTEGER_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC patameter %s REAL %e --\n",yylineno-1,$1,double($3));
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	parameter.set_name($1);
	parameter.set_real($3);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}

	   | STRING_VALUE   STRING   STRING_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC patameter %s STRING %s --\n",yylineno-1,$1,$3);
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	parameter.set_name($1);
	parameter.set_string($3);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}
	   | STRING_VALUE   ENUM  '<' string_values '>'  STRING_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC parameter %s ENUM with default value %s --\n",yylineno-1,$1,$6);
#endif
	//change parameter name to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }

	//change enum string to lower case
	for(unsigned int c=0;c<strlen($6);c++)
        {
          if(isupper($6[c]))
            $6[c]=tolower($6[c]);
        }

	parameter.set_name($1);
	// default ENUM string does not fit ENUM range? user makes mistake with pattern file
	if( parameter.set_enum($6) ) yyerror(dummy, $6);
	card._parameter_map.insert(std::pair<const std::string, const Parser::Parameter>($1,parameter));
	parameter.clear();
}
	   ;

string_values : STRING_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC enum %s --\n",yylineno-1,$1);
#endif
	//change enum string to lower case
	for(unsigned int c=0;c<strlen($1);c++)
        {
          if(isupper($1[c]))
            $1[c]=tolower($1[c]);
        }
	parameter.add_string_pattern($1);
}
              | string_values STRING_VALUE
{
#ifdef VERBOSE
	printf("line %d: -- YACC enum %s --\n",yylineno-1,$2);
#endif
	//change enum string to lower case
	for(unsigned int c=0;c<strlen($2);c++)
        {
          if(isupper($2[c]))
            $2[c]=tolower($2[c]);
        }
	parameter.add_string_pattern($2);
}
	      ;

%%


static int yyerror(void * dummy, const char *s)
{
   std::cerr << "\nYACC report: Pattern file line " << yylineno << " unrecognized word(s): " << s << std::endl;
   return 1;
}

