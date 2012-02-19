%{
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <ctype.h>
#include "input.tab.h"
#include <string>


#define YY_DECL int yylex( void * dummy )

#define KEY_DEFINE      1
#define BOOL_DEFINE     3
#define INT_DEFINE      4
#define REAL_DEFINE     5
#define STRING_DEFINE   6
#define DEFINE_RHS      7

// we have to use this to find parameter type
static std::string keyword;
static int state;
static int error_report(const char *s);
//#define VERBOSE



%}


comment1                      [\/][\/].*
comment2                      [\/][\*]([^\*]|[\*]+[^\*\/])*[\*]*[\/]?
name                          [a-zA-Z\_][a-zA-Z\_0-9\.]*
separator                     [ \r\t\v\f\n]
blank                         [ \t]*
cstring                       [\"]([^\"\\\n]|([\\](.|[\n])))*[\"]
fstring                       [\']([^\'\\\n]|([\\](.|[\n])))*[\']
lineinfo                      \!{cstring}
D       [0-9]
E       [EeDd][+-]?{D}+
SIGN    ([+-]?)

%option noyywrap

%%


^\#[^\n]*                         {/* comment line */      return COMMENT;}

\#[^\n]*                          {/* inline comment */    return COMMENT;}

{comment1}                        {/* c++ style comment */ return COMMENT;}

{lineinfo}                        { yytext[yyleng-1]='\0'; strncpy(yylval.sval, &yytext[2], 1023); return LINEINFO;  }

\n                                {/* new command */ state = 0;}

[ \t\r]*                          {/* use to drop blankspace/tab/cr */}



Title                             {  strncpy(yylval.sval, yytext,1023); return TITLE;}

TRUE         |
True         |
true         |
ON           |
On           |
on                                { yylval.bval = true;    return BOOL_VALUE;}



FALSE         |
False         |
false         |
OFF           |
Off           |
off                                { yylval.bval = false;  return BOOL_VALUE;}


bool                               { state = BOOL_DEFINE;   return BOOL_TYPE;}
int                                { state = INT_DEFINE;    return INT_TYPE;}
double                             { state = REAL_DEFINE;   return REAL_TYPE;}
string                             { state = STRING_DEFINE; return STRING_TYPE;}




^[ \t]*{name}                      {

        // at first, we should drop blank chars
        char * cp = yylval.sval;
        for(unsigned int c=0;c<strlen(yytext)+1;c++)
        {
          if (yytext[c] == '\t' || isspace(yytext[c]))
            continue;
          *cp++ = yytext[c];
        }

        Parser::InputParser * p = (Parser::InputParser *) dummy;
        //does this string match any ID?
        if( state==0 )
        {
          std::string name = yylval.sval;
          if(p->bool_var.count(name))        { return BOOL_ID;       }
          if(p->int_var.count(name))         { return INT_ID;        }
          if(p->real_var.count(name))        { return REAL_ID;       }
          if(p->string_var.count(name))      { return STRING_ID;     }
        }

        //At last, this should be a keyword
        state = KEY_DEFINE;

        //change it to upper case
        for(unsigned int c=0;c<strlen(yylval.sval);c++)
        {
          if(islower(yylval.sval[c]))
            yylval.sval[c]=toupper(yylval.sval[c]);
        }

        // can we find this token in card pattern?
        keyword = yylval.sval;
        if(p->_pattern.ckeck_card(keyword))
        {
          error_report(yytext);
        }
        // return exact card name
        strcpy(yylval.sval, keyword.c_str());
#ifdef VERBOSE
        printf("KEYWORD %s\n",yylval.sval);
#endif
        return KEYWORD;

        //we should never reach here
        error_report(yytext);
}



{name}/{blank}\=                                {

        // A string followed by '='

        //for the state of BOOL_DEFINE, this token should be BOOL_ID
        if (state == BOOL_DEFINE)
        {
          strncpy(yylval.sval, yytext,1023);
#ifdef VERBOSE
        printf("BOOL_ID %s\n",yytext);
#endif
          return BOOL_ID;
        }
        if (state == INT_DEFINE)
        {
          strncpy(yylval.sval, yytext,1023);
#ifdef VERBOSE
        printf("INT_ID %s\n",yytext);
#endif
          return INT_ID;
        }
        if (state == REAL_DEFINE)
        {
          strncpy(yylval.sval, yytext,1023);
#ifdef VERBOSE
        printf("REAL_ID %s\n",yytext);
#endif
          return REAL_ID;
        }
        if (state == STRING_DEFINE)
        {
          strncpy(yylval.sval, yytext,1023);
#ifdef VERBOSE
        printf("STRING_ID %s\n",yytext);
#endif
          return STRING_ID;
        }

        //for the state of KEY DEFINE, it should be parameter
        Parser::InputParser * p = (Parser::InputParser *) dummy;
        if (state == KEY_DEFINE)
        {
          strncpy(yylval.sval, yytext,1023);
          for(unsigned int c=0;c<strlen(yylval.sval);c++)
          {
            if(isupper(yylval.sval[c]))
              yylval.sval[c]=tolower(yylval.sval[c]);
          }

          std::string pname = yylval.sval;
          Parser::ElemType type = p->_pattern.check_parameter_type(keyword, pname);
          strcpy(yylval.sval, pname.c_str());
#ifdef VERBOSE
          printf("PARAMETER_ID %s type %d\n",yylval.sval, type);
#endif
          switch ( type )
          {
            case        Parser::BOOL      :  return BOOL_PARAMETER;
            case        Parser::INTEGER   :  return INT_PARAMETER;
            case        Parser::REAL      :  ;return REAL_PARAMETER;
            case        Parser::STRING    :  ;return STRING_PARAMETER;
            case        Parser::ENUM      :  ;return STRING_PARAMETER;
            default:    break;
          }
          // if user mistakenly input an wrong parameter, however it matches ID
          // YACC will pick this mistake later
        }

        // for a general assign statement, it maybe an ID
        // find if this token belongs to any existing ID
        strncpy(yylval.sval, yytext,1023);
        std::string name = yylval.sval;

        if(p->bool_var.count(name))
        {
#ifdef VERBOSE
              printf("BOOL_ID %s\n",yytext);
#endif
          return BOOL_ID;
        }

        if(p->int_var.count(name))
        {
#ifdef VERBOSE
            printf("INT_ID %s\n",yytext);
#endif
          return INT_ID;
        }
        if(p->real_var.count(name))
        {
#ifdef VERBOSE
            printf("REAL_ID %s\n",yytext);
#endif
          return REAL_ID;
        }
        if(p->string_var.count(name))
        {
#ifdef VERBOSE
            printf("STRING_ID %s\n",yytext);
#endif
          return STRING_ID;
        }

        //we should never reach here
        error_report(yytext);
}


\=                                     { return '=';  }

\+                                     { return PLUS; }

\-                                     { return MINUS; }

\*                                     { return ASTERISK; }

\/                                     { return DIVIDE; }

\(                                     { return '('; }

\)                                     { return ')'; }

\{                                     { return '{'; }

\}                                     { return '}'; }

\[                                     { return '['; }

\]                                     { return ']'; }

\,                                     { return ','; }

\!                                     { return NOT;}

\&\&                                   { return AND;}

\|\|                                   { return OR;}

\>                                     { return GT;}

\>\=                                   { return GE;}

\<                                     { return LT;}

\<=                                    { return LE;}

\=\=                                   { return EQ;}

\!\=                                   { return NE;}

{D}+                                   {
    // NOTE: input deck some times have big number, i.e. 2397610736220310
    // which will be parsed as int. However, int will be overflow for such a big value
    int    i;
    double d;
    sscanf(yytext, "%d", &i);
    sscanf(yytext, "%lf", &d);
    if( i == floor(d+0.5) )
    {
      yylval.ival = i;
      return INT_VALUE;
    }

    yylval.dval = d;
    return REAL_VALUE;
}


{D}+"."{D}*({E})?         |
{D}*"."{D}+({E})?         |
{D}+({E})                       {sscanf(yytext, "%lf", &yylval.dval);return REAL_VALUE;}



{name}                                {

        // find if this token belongs to any exist ID
        strncpy(yylval.sval, yytext,1023);
        std::string name = yylval.sval;

        Parser::InputParser * p = (Parser::InputParser *) dummy;
        if(p->bool_var.count(name))
        {
          return BOOL_ID;
        }

        if(p->int_var.count(name))
        {
          return INT_ID;
        }
        if(p->real_var.count(name))
        {
          return REAL_ID;
        }
        if(p->string_var.count(name))
        {
          return STRING_ID;
        }
#ifdef VERBOSE
        printf("STRING_VALUE %s\n",yylval.sval);
#endif
        // at larst, return as STRING_VALUE;
        return STRING_VALUE;

        //we should never reach here
        error_report(yytext);
}


{cstring}                                    {
        // c style string
        // the only poissable is STRING_VALUE
        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[1],1023);
#ifdef VERBOSE
        printf("STRING_VALUE %s\n",yylval.sval);
#endif
        return STRING_VALUE;
}


{fstring}                                    {
        // fortran style string
        // the only poissable is STRING_VALUE
    yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[1],1023);
#ifdef VERBOSE
        printf("STRING_VALUE %s\n",yylval.sval);
#endif
        return STRING_VALUE;
}



bool\<{name}\>/{blank}\=                     {
       // user defined parameter name
       // the parameter type is bool

        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[strlen("bool<")],1023);
        return UD_BOOL_PARAMETER;
}


int\<{name}\>/{blank}\=                     {
       // user defined parameter name
       // the parameter type is int

        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[strlen("int<")],1023);
        return UD_INT_PARAMETER;
}


double\<{name}\>/{blank}\=                  {
       // user defined parameter name
       // the parameter type is real

        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[strlen("double<")],1023);
        return UD_REAL_PARAMETER;
}


real\<{name}\>/{blank}\=                    {
       // user defined parameter name
       // the parameter type is real

        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[strlen("real<")],1023);
        return UD_REAL_PARAMETER;
}


string\<{name}\>                     {
       // user defined parameter name
       // the parameter type is string

        yytext[yyleng-1]='\0';
        strncpy(yylval.sval, &yytext[strlen("string<")],1023);
        return UD_STRING_PARAMETER;
}


\;        {state = 0; return ';';}


<<EOF>>   {return 0;}

%%

static int error_report(const char * s)
{
  std::cerr  << "ERROR: token "<<s<<" can't be matched."<<std::endl;
  //std::cerr  << "ERROR at " << __FILE__ << ", line " << __LINE__ << std::endl;
  return 1;
}

