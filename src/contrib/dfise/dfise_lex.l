%{

#ifndef YYTOKENTYPE
#include "dfise_parser.tab.h"
#endif

/*#define VERBOSE*/

%}

D   [0-9]
E   [Ee][+-]?{D}+
L   [a-zA-Z_]
SIGN    ([+-]?)


%%

DF\-ISE                            {
#ifdef VERBOSE
    ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return DFISE;
}


text                     |
binary                             {
#ifdef VERBOSE
    printf("FILE_FORMAT:");ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return FILE_FORMAT;
}


Info         |
Data         |
Dataset      |
CoordSystem  |
Vertices     |
Edges        |
Faces        |
Locations    |
Elements     |
Region       |
Values                            {
#ifdef VERBOSE
    printf("KEYWORD:");ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return KEYWORD;
}


version      |
type         |
dimension    |
nb\_vertices |
nb\_edges    |
nb\_faces    |
nb\_elements |
nb\_regions  |
regions      |
materials    |
datasets     |
functions                        {
#ifdef VERBOSE
    printf("PARAMETER:");ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return KEYWORD;
}



translate    |
transform    |
material                         {
#ifdef VERBOSE
    printf("PARAMETER:");ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return KEYWORD;
}


function     |
location     |
validity                          {
#ifdef VERBOSE
    printf("PARAMETER:");ECHO;
#endif
    strncpy(yylval.sval,yytext,255);
    return KEYWORD;
}


{SIGN}{D}+                         {
#ifdef VERBOSE
    printf("INTEGER:");ECHO;
#endif
    // NOTE: DF-ISE some times have float number i.e 2397610736220310
    // which will be parsed as int. However, int will be overflow for such a big value
    int    i;
    double d;
    sscanf(yytext, "%d", &i);
    sscanf(yytext, "%lf", &d);
    if( i == floor(d+0.5) )
    {
      yylval.ival = i;
      return INTEGER;
    }

    yylval.dval = d;
    return FLOAT;
}


{SIGN}{D}+"."{D}*({E})? |
{SIGN}{D}*"."{D}+({E})? |
{SIGN}{D}+({E})                    {
#ifdef VERBOSE
     printf("FLOAT:");ECHO;
#endif
     sscanf(yytext, "%lf", &yylval.dval);
     return FLOAT;
}



[a-zA-Z0-9\_][a-zA-Z0-9\.\_\+\&]*  {
#ifdef VERBOSE
        printf("STRING:%s",yytext);
#endif
        strncpy(yylval.sval,yytext,255);
        return STRING;
}



\"[a-zA-Z0-9\.\_\ \-\+\&]*\"       {
        yytext[yyleng-1]='\0';
#ifdef VERBOSE
        printf("STRING:%s",&yytext[1]);
#endif
        strncpy(yylval.sval, &yytext[1],255);
        return STRING;
}



\{               { return '{'; }
\}               { return '}'; }
\[               { return '['; }
\]               { return ']'; }
\(               { return '('; }
\)               { return ')'; }
\=               { return '='; }

\#[^\n]*                     {
/* comment */
#ifdef VERBOSE
    ECHO;
#endif
}

[ \t\r\n]                      {
/* use to drop blankspace/tab */
#ifdef VERBOSE
    ECHO;
#endif
}



%%

