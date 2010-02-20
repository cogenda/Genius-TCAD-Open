%{
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <ctype.h>
#include "pattern.tab.h"

//#define VERBOSE


%}
comment1                      [\/][\/].*
comment2                      [\/][\*]([^\*]|[\*]+[^\*\/])*[\*]*[\/]?
comment3                      ^\#[^\n]*
comment4                      \#[^\n]*

D       [0-9]
E       [EeDd][+-]?{D}+
SIGN    ([+-]?)

%option noyywrap

%%

{comment1}                 |
{comment2}                 |
{comment3}                 |
{comment4}                          {
/* comment line */
#ifdef VERBOSE
        printf("COMMENT LINE:");ECHO;
#endif
}



BOOL                     |
LOG                               {
#ifdef VERBOSE
        ECHO;
#endif
        return BOOL;
}


INTEGER                  |
INT                               {
#ifdef VERBOSE
        ECHO;
#endif
        return INTEGER;
}

REAL                     |
NUM                               {
#ifdef VERBOSE
        ECHO;
#endif
        return REAL;
}

STRING                            {
#ifdef VERBOSE
        ECHO;
#endif
        return STRING;
}

ENUM                              {
#ifdef VERBOSE
        ECHO;
#endif
        return ENUM;
}

{SIGN}{D}+                        {
        sscanf(yytext, "%d", &yylval.ival);
#ifdef VERBOSE
        printf("int_value: %d",yylval.ival);
#endif
        return INTEGER_VALUE;
}


{SIGN}{D}+"."{D}*({E})?     |
{SIGN}{D}*"."{D}+({E})?     |
{SIGN}{D}+({E})                   {
        sscanf(yytext, "%lf", &yylval.dval);
#ifdef VERBOSE
        printf("double_value: %lf",yylval.dval);
#endif
        return REAL_VALUE;
}



TRUE         |
True         |
true         |
ON           |
On           |
on           |
T            |
t                                 {
        yylval.bval = true;
#ifdef VERBOSE
        printf("BOOL True:%s",yytext);
#endif
        return BOOL_VALUE;
}



FALSE         |
False         |
false         |
OFF           |
Off           |
off           |
F             |
f                                 {
        yylval.bval = false;
#ifdef VERBOSE
        printf("BOOL False:%s",yytext);
#endif
        return BOOL_VALUE;
}


[a-zA-Z0-9\_][a-zA-Z0-9\.\_]*     {
#ifdef VERBOSE
        printf("STRING:%s",yytext);
#endif
        strncpy(yylval.sval,yytext,31);
        return STRING_VALUE;
}



[\"]([^\"\\\n]|([\\](.|[\n])))*[\"]        {
        yytext[yyleng-1]='\0';
#ifdef VERBOSE
        printf("STRING:%s",&yytext[1]);
#endif
        strncpy(yylval.sval, &yytext[1],31);
        return STRING_VALUE;
}


\{                                {
#ifdef VERBOSE
        printf("{");
#endif
        return '{';
}

\}                                {
#ifdef VERBOSE
        printf("}");
#endif
        return '}';
}

\<                                {
#ifdef VERBOSE
        printf("<");
#endif
        return '<';
}

\>                                {
#ifdef VERBOSE
        printf(">");
#endif
        return '>';
}

[ \t\r\n]*                          {
/* use to drop blankspace/tab/cr */
#ifdef VERBOSE
        ECHO;
#endif

}



%%

