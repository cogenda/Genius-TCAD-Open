%{
#include <stdio.h>
#include <string.h>


extern int yylineno;
extern int yylex();
static int yyerror(void * dummy, const char *s);



static struct Parser::Card card;
static Parser::Parameter   parameter;

/*#define VERBOSE*/

%}

%start inputfile

%union{
    bool   bval;
    int    ival;
    double dval;
    char   sval[1024];
}

%parse-param { void * dummy}
%lex-param { void * dummy}

%token <sval> LINEINFO TITLE KEYWORD  BOOL_PARAMETER REAL_PARAMETER INT_PARAMETER STRING_PARAMETER
%token <sval> UD_BOOL_PARAMETER UD_REAL_PARAMETER UD_INT_PARAMETER UD_STRING_PARAMETER
%token INT_TYPE REAL_TYPE BOOL_TYPE STRING_TYPE
%token <sval> INT_ID REAL_ID BOOL_ID STRING_ID
%token <ival> INT_VALUE
%token <dval> REAL_VALUE
%token <sval> STRING_VALUE
%token <bval> BOOL_VALUE
%token ENDFILE

%token  BLANK COMMENT BAD_WORD
%token  GE GT LE LT EQ NE
%token  NOT OR AND
%token  PLUS MINUS ASTERISK DIVIDE MOD


%nonassoc LOWEST
%nonassoc ELSE
%left OR
%left AND
%left EQ NE
%left LE LT GE GT
%left PLUS MINUS
%left MOD ASTERISK DIVIDE
%nonassoc UNARYLOW
%nonassoc UNARYHIGH
%nonassoc LEFTSQUARE LEFTPAREN
%nonassoc HIGHEST

%type  <bval> bool_expr
%type  <ival> int_expr
%type  <dval> real_expr
%type  <sval> string_expr

%%
inputfile
    :   TITLE STRING_VALUE
    |   command
    |   inputfile command
    ;

command
     :   statement

     |   BOOL_TYPE     BOOL_ID    '='  bool_expr    ';'
         {
            Parser::InputParser * p = (Parser::InputParser *) dummy;
            if(p->bool_var.count($2))
            {
              std::cerr<<"Error: bool value "<<$2<<" has already exist."<<std::endl;
              yyerror(dummy,$2);
            }
            p->bool_var.insert(std::pair< const std::string, bool >($2,$4));
         }
     |   INT_TYPE      INT_ID     '='  int_expr     ';'
         {
            Parser::InputParser * p = (Parser::InputParser *) dummy;
            if(p->int_var.count($2))
            {
              std::cerr<<"Error: int value "<<$2<<" has already exist."<<std::endl;
              yyerror(dummy,$2);
            }
            p->int_var.insert(std::pair< const std::string, int >($2,$4));
         }
     |   REAL_TYPE     REAL_ID    '='  real_expr    ';'
         {
            Parser::InputParser * p = (Parser::InputParser *) dummy;
            if(p->real_var.count($2))
            {
              std::cerr<<"Error: real value "<<$2<<" has already exist."<<std::endl;
              yyerror(dummy,$2);
            }
            p->real_var.insert(std::pair< const std::string, double >($2,$4));
         }
     |   STRING_TYPE   STRING_ID  '='  string_expr  ';'
         {
            Parser::InputParser * p = (Parser::InputParser *) dummy;
            if(p->string_var.count($2))
            {
              std::cerr<<"Error: string value "<<$2<<" has already exist."<<std::endl;
              yyerror(dummy,$2);
            }
            p->string_var.insert(std::pair< const std::string, std::string >($2,$4));
         }

    |   LINEINFO
        {
           //LINEINFO is always at the end of card
           if( !card.key().empty() )
           {
             card.set_lineno(yylineno);
             card.set_fileline($1);
             Parser::InputParser * p = (Parser::InputParser *) dummy;
             p->_card_list.push_back(card);
             p->_card_map.insert(std::pair<const std::string, Parser::Card>(card.key(), card));
             card.clear();
           }
        }

    |   COMMENT

    |   BLANK
    ;


statement : KEYWORD
{
#ifdef VERBOSE
        printf("line %d: -- YACC command %s --\n",yylineno,$1);
#endif
           card.set_key($1);
}
          | KEYWORD parameters
{
#ifdef VERBOSE
        printf("line %d: -- YACC command %s --\n",yylineno,$1);
#endif
           card.set_key($1);
}
          |   expr ';'
          ;

expr      :
          | BOOL_ID   '=' bool_expr
          | INT_ID    '=' int_expr
          | REAL_ID   '=' real_expr
          | STRING_ID '=' string_expr
          ;



bool_expr : BOOL_VALUE                             {$$=$1;}
          | BOOL_ID
            {
              Parser::InputParser * p = (Parser::InputParser *) dummy;
              if(!p->bool_var.count($1))
              {
                   std::cerr  << "ERROR: bool ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
              }
              $$=(*(p->bool_var.find($1))).second;
            }
          | bool_expr AND bool_expr                {$$=($1 && $3);}
          | bool_expr OR bool_expr                 {$$=($1 || $3);}
          | bool_expr EQ bool_expr                 {$$=($1 == $3);}
          | bool_expr NE bool_expr                 {$$=($1 != $3);}
          | '('  bool_expr  ')'                    {$$=$2;}
          | NOT bool_expr %prec UNARYLOW           {$$=!$2;}
          | real_expr LT real_expr                 {$$=($1 < $3);}
          | real_expr GT real_expr                 {$$=($1 > $3);}
          | real_expr LE real_expr                 {$$=($1 <= $3);}
          | real_expr GE real_expr                 {$$=($1 >= $3);}
          | real_expr EQ real_expr                 {$$=($1 == $3);}
          | real_expr NE real_expr                 {$$=($1 != $3);}
          ;

int_expr  : INT_VALUE                              {$$=$1;}
          | INT_ID
            {
              Parser::InputParser * p = (Parser::InputParser *) dummy;
              if(!p->int_var.count($1))
              {
                   std::cerr  << "ERROR: int ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
              }
              $$=(*(p->int_var.find($1))).second;
            }
          | '('  int_expr  ')'                     {$$=$2;}
          | int_expr      PLUS     int_expr        {$$=$1+$3;}
          | int_expr      MINUS    int_expr        {$$=$1-$3;}
          | int_expr    ASTERISK   int_expr        {$$=$1*$3;}
          | int_expr     DIVIDE    int_expr        {$$=$1/$3;}
          | MINUS int_expr  %prec UNARYLOW         {$$=-$2;}
          | PLUS  int_expr  %prec UNARYLOW         {$$= $2;}
          ;

real_expr: REAL_VALUE                              {$$=$1;}
         | INT_VALUE                               {$$=$1;}
         | REAL_ID
           {
              Parser::InputParser * p = (Parser::InputParser *) dummy;
              if(!p->real_var.count($1))
              {
                   std::cerr  << "ERROR: real ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
              }
              $$=(*(p->real_var.find($1))).second;
           }
         | INT_ID
           {
              Parser::InputParser * p = (Parser::InputParser *) dummy;
              if(!p->int_var.count($1))
              {
                   std::cerr  << "ERROR: int ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
              }
              $$=(*(p->int_var.find($1))).second;
            }
         | '(' real_expr ')'                       {$$=$2;}
         | real_expr      PLUS     real_expr       {$$=$1+$3;}
         | real_expr      MINUS    real_expr       {$$=$1-$3;}
         | real_expr    ASTERISK  real_expr        {$$=$1*$3;}
         | real_expr     DIVIDE   real_expr        {$$=$1/$3;}
         | MINUS real_expr  %prec UNARYLOW         {$$=-$2;}
         | PLUS  real_expr  %prec UNARYLOW         {$$= $2;}
         ;

string_expr :  STRING_VALUE
               {
                 strncpy($$,$1,1023);
               }
            |  STRING_ID
               {
                 Parser::InputParser * p = (Parser::InputParser *) dummy;
                 if(!p->string_var.count($1))
                 {
                   std::cerr  << "ERROR: string ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
                 }
                 std::string s = (*(p->string_var.find($1))).second;
                 strncpy($$,s.c_str(),1023);
               }
            |  INT_ID
               {
                 Parser::InputParser * p = (Parser::InputParser *) dummy;
                 if(!p->int_var.count($1))
                 {
                   std::cerr  << "ERROR: int ID "<<$1<<" can't be find in variable table."<<std::endl;
                   yyerror(dummy, $1);
                 }
                 int n=(*(p->int_var.find($1))).second;
                 sprintf($$,"%d",n);
               }
            |  string_expr PLUS string_expr
               {
                  std::string s1 = $1;
                  s1 = s1+ $3;
                  strncpy($$,s1.c_str(),1023);
               }
            ;

parameters  :  parameter
            |  parameters parameter
            ;

parameter
     : REAL_PARAMETER '=' real_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC patameter %s REAL %e --\n",yylineno,$1,$3);
#endif
        parameter.set_name($1);
        parameter.set_real($3);
        card.insert(parameter);
        parameter.clear();
}
     | UD_REAL_PARAMETER '=' real_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC user defined patameter %s REAL %e --\n",yylineno,$1,$3);
#endif
        parameter.set_name($1);
        parameter.set_real($3);
        parameter.do_not_check_me = true;
        card.insert(parameter);
        parameter.clear();
}

     | INT_PARAMETER '=' int_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC patameter %s INT %d --\n",yylineno,$1,int($3));
#endif
        parameter.set_name($1);
        parameter.set_int($3);
        card.insert(parameter);
        parameter.clear();
}
     | UD_INT_PARAMETER '=' int_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC user defined patameter %s INT %d --\n",yylineno,$1,int($3));
#endif
        parameter.set_name($1);
        parameter.set_int($3);
        parameter.do_not_check_me = true;
        card.insert(parameter);
        parameter.clear();
}

     | STRING_PARAMETER '=' string_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC patameter %s STRING %s --\n",yylineno,$1,$3);
#endif
        parameter.set_name($1);
        parameter.set_string($3);
        card.insert(parameter);
        parameter.clear();
}
     | UD_STRING_PARAMETER '=' string_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC user defined patameter %s STRING %s --\n",yylineno,$1,$3);
#endif
        parameter.set_name($1);
        parameter.set_string($3);
        parameter.do_not_check_me = true;
        card.insert(parameter);
        parameter.clear();
}

     | STRING_PARAMETER
{
#ifdef VERBOSE
        printf("line %d: -- YACC patameter %s STRING --\n",yylineno,$1);
#endif
        parameter.set_name($1);
        parameter.set_string("");
        card.insert(parameter);
        parameter.clear();
}
     | BOOL_PARAMETER '=' bool_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC patameter %s BOOL %d --\n",yylineno,$1,int($3));
#endif
        parameter.set_name($1);
        parameter.set_bool($3);
        card.insert(parameter);
        parameter.clear();
}
     | UD_BOOL_PARAMETER '=' bool_expr
{
#ifdef VERBOSE
        printf("line %d: -- YACC user defined patameter %s BOOL %d --\n",yylineno,$1,int($3));
#endif
        parameter.set_name($1);
        parameter.set_bool($3);
        parameter.do_not_check_me = true;
        card.insert(parameter);
        parameter.clear();
}
     ;

%%

static int yyerror(void * dummy, const char *)
{
   std::cerr << "\nYACC report: Input file line " << yylineno << " has unrecognized word(s)." << std::endl;
   std::cerr << "Please check the syntax and try again." << std::endl;
   return 1;
}

