%{
#include "tif_parser.tab.h"
/*#define VERBOSE*/

%}

D       [0-9]
E       [Ee][+-]?{D}+
L       [a-zA-Z_]
SIGN    ([+-]?)

RC      ^[ \t]*c[ \t]+
RE      ^[ \t]*e[ \t]+
RR      ^[ \t]*r[ \t]+
RB      ^[ \t]*b[ \t]+
RI      ^[ \t]*i[ \t]+
RJ      ^[ \t]*j[ \t]+
RT      ^[ \t]*t[ \t]+
RQ      ^[ \t]*q[ \t]+
RS      ^[ \t]*s[ \t]+
RN      ^[ \t]*n[ \t]+
COMP    ^[ \t]*comp[ \t]+
MAT     ^[ \t]*mo[ \t]+blank[ \t]+mat[ \t]+

%option noyywrap
%Start RECORD MATERIAL IGNORE
%%

\#.*                       {
#ifdef VERBOSE
        printf("COMMENT:");ECHO;
#endif
}


{RC}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Coordinate:");ECHO;
#endif
        return Coordinate;
}


{RE}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Edge:");ECHO;
#endif
        return Edge;
}


{RR}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Region:");ECHO;
#endif
        return Region;
}


{RB}        {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Boundary:");ECHO;
#endif
        return Boundary;
}


{RI}        {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Interface:");ECHO;
#endif
        return Interface;
}


{RJ}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD InterfaceBD:");ECHO;
#endif
        return InterfaceBD;
}


{RT}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Triangle:");ECHO;
#endif
        return Triangle;
}

{RQ}         {
  BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Quadrangle:");ECHO;
#endif
        return Quadrangle;
}


{RS}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Solution:");ECHO;
#endif
        return Solution;
}


{RN}         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD NodalSolution:");ECHO;
#endif
        return NodalSolution;
}

{COMP}       {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("RECORD Component:");ECHO;
#endif
        return Component;
}

<RECORD>{SIGN}{D}+                 {
#ifdef VERBOSE
        printf("interger");ECHO;
#endif
        sscanf(yytext, "%d", &yylval.ival);
        return Integer;
}

<RECORD>{SIGN}{D}+"."{D}*({E})? |
<RECORD>{SIGN}{D}*"."{D}+({E})? |
<RECORD>{SIGN}{D}+({E})           {
#ifdef VERBOSE
         printf("float:");ECHO;
#endif
         sscanf(yytext, "%lf", &yylval.dval);
         return Float;
}


<RECORD>[a-zA-Z\_]+[a-zA-Z0-9\(\)\.\/\-\_]*     {
#ifdef VERBOSE
        printf("string:"); ECHO;
#endif
        strcpy(yylval.sval, yytext);
        return String;
}

\n                           {
BEGIN 0;
#ifdef VERBOSE
        ECHO;
#endif
}

[ \t\r]                      {
/* use to drop blankspace/tab */
#ifdef VERBOSE
        ECHO;
#endif
}


^[ \t]*ds[ \t]+         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("Data string:");ECHO;
#endif
        return DataString;
}

^[ \t]*di[ \t]+         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("Data string:");ECHO;
#endif
        return Dataname;
}

\+[ \t]*\n              {
#ifdef VERBOSE
        printf("continue line:");ECHO;
#endif
        return CLine;
}

^[ \t]*du[ \t]+         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("Data Unin:");ECHO;
#endif
        return DataUnit;
}

^[ \t]*dd[ \t]+         {
BEGIN  RECORD;
#ifdef VERBOSE
        printf("Data DataArray:");ECHO;
#endif
        return DataArray;
}

{MAT} {BEGIN MATERIAL; return MAT;}

<MATERIAL>{SIGN}{D}+[ \t]+                {
        sscanf(yytext, "%d", &yylval.ival);
        return Integer;
}
<MATERIAL>\= { return '='; }
<MATERIAL>{SIGN}{D}+"."{D}*({E})? |
<MATERIAL>{SIGN}{D}*"."{D}+({E})? |
<MATERIAL>{SIGN}{D}+({E})                 {
         sscanf(yytext, "%lf", &yylval.dval);
         return Float;
}
<MATERIAL>[a-zA-Z0-9\(\)\.\/\-]+          {
        strcpy(yylval.sval, yytext);
        return PARAMETER;
}


^[ \t]*[a-zA-Z0-9]+    {
/* unrecognized chars */
BEGIN IGNORE;
#ifdef VERBOSE
      printf("enter unwanted ");ECHO;
#endif
}

<IGNORE>.* {
#ifdef VERBOSE
        ECHO;
        printf("line %d \n",yylineno);
#endif
}

%%

