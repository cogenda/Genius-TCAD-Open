%{

extern int yylineno;
extern int yylex();
int yyerror(BLOCK *, char *s);

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
    std::vector<TOKEN *> * tokens;
    TOKEN * token;
   }

%token <sval> DFISE   FILE_FORMAT
%token <ival> INTEGER
%token <dval> FLOAT
%token <sval> STRING KEYWORD PARAMETER

%type <bval> block body bodyitem
%type <token> value
%type <tokens> values data

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
            $$->add_values(*$1);
            delete $1;
}
         | KEYWORD '=' data
{
           $$ = new BLOCK;
           $$->add_parameter($1, *$3);
           delete $3;
}
         | KEYWORD '=' KEYWORD
{
           $$ = new BLOCK;
           TOKEN * new_token = new TOKEN;
           new_token->token_type = TOKEN::string_token;
           new_token->value = new std::string;
           *((std::string*)new_token->value) = $3;
           std::vector<TOKEN *> tokens;
           tokens.push_back(new_token);
           $$->add_parameter($1, tokens);
}
         ;



data     : value
{
         $$ = new std::vector<TOKEN *>;
         $$->push_back($1);
}
         | '[' values ']'
{
         $$ = $2;
}
         ;

values   :  value
{
         $$ = new std::vector<TOKEN *>;
         $$->push_back($1);
}
         |  values value
{
         $$ = $1;
         $$->push_back($2);
}
         ;


value    : INTEGER
{
         // read integer data, insert into vector TOKENS
         $$ = new TOKEN;
         $$->token_type = TOKEN::int_token;
         $$->value = new int;
         *((int*)$$->value) = $1;
}
         | FLOAT
{
         // read float data, insert into vector TOKENS
         $$ = new TOKEN;
         $$->token_type = TOKEN::float_token;
         $$->value = new double;
         *((double*)$$->value) = $1;
}
         | STRING
{
         // read string data, insert into vector TOKENS
         $$ = new TOKEN;
         $$->token_type = TOKEN::string_token;
         $$->value = new std::string;
         *((std::string*)$$->value) = $1;
}
         ;

%%

int yyerror(BLOCK *, char *)
{
   printf("\nline %d unrecognized chars %s \n",yylineno, yylval.sval);
   return 0;
}



