%{

extern int yylineno;
extern int yylex();
int yyerror(BLOCK *, char *s);

std::vector<TOKEN *> TOKENS;


//#define VERBOSE

%}
%start file

%parse-param {BLOCK * ise_block}

%union  {
    int    ival;
    double dval;
    char   cval;
    char   sval[256];
    BLOCK * bval;
   }

%token <sval> DFISE   FILE_FORMAT
%token <ival> INTEGER
%token <dval> FLOAT
%token <sval> STRING KEYWORD PARAMETER

%type <bval> block body bodyitem

%%

file     : DFISE FILE_FORMAT blocks
{
           ise_block->set_keyword($1);
}
         ;


blocks   : block
{
           /* top level block */
           ise_block->add_sub_block($1);
}
         | blocks block
{
           /* top level block */
           ise_block->add_sub_block($2);
}
         ;


block    : KEYWORD '{' body '}'
{
#ifdef VERBOSE
    printf("block1:%s {body}\n", $1);
#endif
         /* new block */
         $$ = new BLOCK($1);
         $$->append(*$3);
         delete $3;
}
         | KEYWORD '(' INTEGER ')' '{' body '}'
{
#ifdef VERBOSE
    printf("block2:%s (%d) {body}\n", $1, $3);
#endif
         /* new block */
         $$ = new BLOCK($1);
         $$->set_index($3);
         $$->append(*$6);
         delete $6;
}
         | KEYWORD '(' STRING ')' '{' body '}'
{
#ifdef VERBOSE
    printf("block3:%s (%s) {body}\n", $1, $3);
#endif
         /* new block */
         $$ = new BLOCK($1);
         $$->set_label($3);
         $$->append(*$6);
         delete $6;
}
         ;

body     :  bodyitem
{
         $$ = new BLOCK;
         $$->append(*$1);
         delete $1;
}
         |  body bodyitem
{
         $$ = $1;
         $$->append(*$2);
         delete $2;
}
         ;

bodyitem :  block
{
            $$ = new BLOCK;
            $$->add_sub_block($1);
}
         |  data
{
            $$ = new BLOCK;
            /* data is pushed into TOKENS vector */
            $$->add_values(TOKENS);
            TOKENS.clear();
}
         | KEYWORD '=' data
{
           $$ = new BLOCK;
           $$->add_parameter($1, TOKENS);
           TOKENS.clear();
}
         | KEYWORD '=' KEYWORD
{
           $$ = new BLOCK;
           TOKEN * new_token = new TOKEN;
           new_token->token_type = TOKEN::string_token;
           new_token->value = new std::string;
           *((std::string*)new_token->value) = $3;
           TOKENS.push_back(new_token);
           $$->add_parameter($1, TOKENS);
           TOKENS.clear();
}
         ;



data     : value
         | '[' values ']'
         ;

values   :  value
         |  values value
         ;

value    : INTEGER
{

         // read integer data, insert into vector TOKENS
         TOKEN * new_token = new TOKEN;
         new_token->token_type = TOKEN::int_token;
         new_token->value = new int;
         *((int*)new_token->value) = $1;
         TOKENS.push_back(new_token);

}
         | FLOAT
{

         // read float data, insert into vector TOKENS
         TOKEN * new_token = new TOKEN;
         new_token->token_type = TOKEN::float_token;
         new_token->value = new double;
         *((double*)new_token->value) = $1;
         TOKENS.push_back(new_token);

}
         | STRING
{

         // read string data, insert into vector TOKENS
         TOKEN * new_token = new TOKEN;
         new_token->token_type = TOKEN::string_token;
         new_token->value = new std::string;
         *((std::string*)new_token->value) = $1;
         TOKENS.push_back(new_token);

}
         ;

%%

int yyerror(BLOCK *, char *)
{
   printf("\nline %d unrecognized chars %s \n",yylineno, yylval.sval);
   return 0;
}


