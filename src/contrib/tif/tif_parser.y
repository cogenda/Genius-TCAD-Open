%{

extern int yylineno;
extern int yylex();
int yyerror(char *s);

/*#define VERBOSE*/

/* global data storage */
std::vector<Node_t>      node_array;
std::vector<Edge_t>      edge_array;
std::vector<Region_t>    region_array;
std::vector<Interface_t> interface_array;
std::vector<Tri_t>       tri_array;
std::vector<Quad_t>      quad_array;
std::vector<Parameter_t> parameter_array;
std::vector<Component_t> component_array;
SolHead_t                sol_head;
std::vector<SolData_t>   sol_data;

TranSol_t                tran_sol;

/* local temporary variable */
static Node_t      tmpnode;
static Edge_t      tmpedge;
static Region_t    tmpregion;
static Interface_t tmpinterface;
static Tri_t       tmptri;
static Quad_t      tmpquad;
static Parameter_t tmpparameter;
static Component_t tmpcomponent;
static SolData_t   tmpdata;
static std::string tmpstring;

%}

%start tiffile

%union  {
    int    ival;
    double dval;
    char   cval;
    char   sval[256];
   }


%token <cval> Coordinate Edge Region Boundary Interface
%token <cval> InterfaceBD Triangle Quadrangle Solution NodalSolution Component
%token <ival> DataString Dataname DataUnit DataArray CLine MAT
%token <ival> Integer
%token <dval> Float
%token <sval> String PARAMETER
%token <ival> Ignore
%type  <dval> float
%type  <ival> integer

%%
tiffile
	:	record
	|	tiffile record
	;

float   : Integer {$$ = double($1);}
	| Float   {$$ = $1;}
	;

integer : Integer {$$ = $1;}
	| Float   {$$ = int($1+0.5);}
	;

solname :  solname String
{
#ifdef VERBOSE
	printf("YACC solname %s\n",$2);
#endif
	tmpstring=$2;
	sol_head.sol_name_array.push_back(tmpstring);
}
	|  String
{
#ifdef VERBOSE
	printf("YACC solname %s\n",$1);
#endif
	tmpstring=$1;
	sol_head.sol_name_array.push_back(tmpstring);
}
	;



value   : value float  {tmpdata.data_array.push_back($2);}
	| float        {tmpdata.data_array.push_back($1);}
	;

parameters : parameters PARAMETER '=' float {
		tmpstring = $2;
		tmpparameter.parameter_name_array.push_back(tmpstring);
		tmpparameter.parameter_value_array.push_back($4);
	     }
	   | PARAMETER '=' float            {
		tmpstring = $1;
		tmpparameter.parameter_name_array.push_back(tmpstring);
		tmpparameter.parameter_value_array.push_back($3);
	     }
	   ;

boundary: boundary Boundary integer
{
#ifdef VERBOSE
	 printf("YACC Boundary %d \n",$3);
#endif
	 tmpregion.boundary.push_back($3-1);

}
	| Boundary integer
{
#ifdef VERBOSE
	 printf("YACC Boundary %d \n",$2);
#endif
	 tmpregion.boundary.push_back($2-1);
}
	;


interfacebd: interfacebd InterfaceBD integer
{
#ifdef VERBOSE
	 printf("YACC InterfaceBD %d \n",$3);
#endif
	 tmpinterface.boundary.push_back($3-1);
}
	 | InterfaceBD integer
{
#ifdef VERBOSE
	 printf("YACC InterfaceBD %d \n",$2);
#endif
	 tmpinterface.boundary.push_back($2-1);
}
	 ;

namestring : namestring String CLine	 { tmpstring=$2; tran_sol.data_name_array.push_back(tmpstring);}
	   | namestring String           { tmpstring=$2; tran_sol.data_name_array.push_back(tmpstring);}
	   | String                      { tmpstring=$1; tran_sol.data_name_array.push_back(tmpstring);}
	   ;

unitstring : unitstring String CLine     { tmpstring=$2; tran_sol.data_unit_array.push_back(tmpstring);}
 	   | unitstring String 		 { tmpstring=$2; tran_sol.data_unit_array.push_back(tmpstring);}
	   | String 			 { tmpstring=$1; tran_sol.data_unit_array.push_back(tmpstring);}
           ;

datalist : datalist float CLine		 { tran_sol.data_value_array.push_back($2);}
	 | datalist float		 { tran_sol.data_value_array.push_back($2);}
	 | float			 { tran_sol.data_value_array.push_back($1);}
	 ;

record	: Coordinate  integer  float  float  float
{
#ifdef VERBOSE
	printf("YACC Coordinate %d %e %e %e\n",$2,$3,$4,$5);
#endif
	tmpnode.index = $2-1;
	tmpnode.x = $3;
	tmpnode.y = $4; /*invert y axis*/
	tmpnode.h = $5;
	node_array.push_back(tmpnode);
}

	| Edge  integer  integer  integer  integer
{
#ifdef VERBOSE
	 printf("YACC edge %d %d %d %d\n",$2,$3,$4,$5);
#endif
	 tmpedge.index  = $2-1;
	 tmpedge.point1 = $3 < $4 ? $3-1 : $4-1 ;
	 tmpedge.point2 = $3 > $4 ? $3-1 : $4-1 ;
	 tmpedge.bcode  = $5;
	 tmpedge.bc_index = -1024;
	 edge_array.push_back(tmpedge);
}

        | Region integer String String
{
#ifdef VERBOSE
    printf("YACC Region %d %s %s\n",$2,$3,$4);
#endif
	 tmpregion.index = $2-1;
	 tmpregion.type  = $3;
	 tmpregion.name  = $4;
	 //sprintf(tmpregion.name,"%s%d", $4, tmpregion.index);
	 region_array.push_back(tmpregion);
	 tmpregion.boundary.clear();
}

	| Region integer String String	boundary
{
#ifdef VERBOSE
	 printf("YACC Region %d %s %s\n",$2,$3,$4);
#endif
	 tmpregion.index = $2-1;
	 tmpregion.type  = $3;
	 tmpregion.name  = $4;
	 //sprintf(tmpregion.name,"%s%d", $4, tmpregion.index);
	 region_array.push_back(tmpregion);
	 tmpregion.boundary.clear();
}

	| Region integer String integer boundary
{
/* some tif file use an integer value instead of region name */
#ifdef VERBOSE
	 printf("YACC Region %d %s %d\n",$2,$3,$4);
#endif
	 tmpregion.index = $2-1;
	 tmpregion.type  = $3;
         char buf[32];
	 sprintf(buf,"r%d", $4);
	 tmpregion.name = buf;
	 region_array.push_back(tmpregion);
	 tmpregion.boundary.clear();
}

	| Interface integer String String integer interfacebd
{
#ifdef VERBOSE
	 printf("YACC Interface %d %s %s region %d\n",$2,$3,$4,$5);
#endif
	 tmpinterface.index  = $2-1;
	 tmpinterface.type   = $3;
	 tmpinterface.name   = $4;
	 tmpinterface.region = $5-1;
	 interface_array.push_back(tmpinterface);
	 tmpinterface.boundary.clear();
}

	| Interface integer String String integer
{
#ifdef VERBOSE
	 printf("YACC Electrode %d %s %s region %d\n",$2,$3,$4,$5);
#endif
	 tmpinterface.index  = $2-1;
	 tmpinterface.type   = $3;
	 tmpinterface.name   = $4;
	 tmpinterface.region = $5-1;
	 interface_array.push_back(tmpinterface);
	 tmpinterface.boundary.clear();
}
        | Interface integer String integer  integer interfacebd
{
/* some tif file use an integer value instead of region name */
#ifdef VERBOSE
	 printf("YACC Electrode %d %s %d region %d\n",$2,$3,$4,$5);
#endif
	 tmpinterface.index = $2-1;
	 tmpinterface.type  = $3;
         char buf[32];
	 sprintf(buf, "r%d", $4);
	 tmpregion.name      = buf;
	 tmpinterface.region = $5-1;
	 interface_array.push_back(tmpinterface);
	 tmpinterface.boundary.clear();
}
	| Interface integer String integer  integer
{
#ifdef VERBOSE
	 printf("YACC Electrode %d %s %d region %d\n",$2,$3,$4,$5);
#endif
	 tmpinterface.index = $2-1;
	 tmpinterface.type  = $3;
	 char buf[32];
	 sprintf(buf, "r%d", $4);
	 tmpregion.name      = buf;
	 tmpinterface.region = $5-1;
	 interface_array.push_back(tmpinterface);
	 tmpinterface.boundary.clear();
}
	| Triangle integer integer integer integer integer integer integer integer
{
#ifdef VERBOSE
	 printf("YACC Triangle node : \t%d\t%d\t%d\n",$4,$5,$6);
#endif
	 tmptri.index = $2-1;
	 tmptri.region = $3-1;
	 tmptri.c1 = $4-1;
	 tmptri.c2 = $5-1;
	 tmptri.c3 = $6-1;
	 tmptri.t1 = $7-1;
	 tmptri.t2 = $8-1;
	 tmptri.t3 = $9-1;
	 if(tmptri.c1==-1 || tmptri.c2==-1 || tmptri.c3==-1)
	 {
	   printf("Warning, degradation triangle meet at index %d. Ignored.\n",tmptri.index);
	 }
         else
           tri_array.push_back(tmptri);
}
         | Quadrangle integer integer integer integer integer integer
{
#ifdef VERBOSE
    printf("YACC Quadrangle node : \t%d\t%d\t%d\n",$4,$5,$6);
#endif
         tmpquad.index = $2-1;
         tmpquad.region = $3-1;
         tmpquad.c1 = $4-1;
         tmpquad.c2 = $5-1;
         tmpquad.c3 = $6-1;
         tmpquad.c4 = $7-1;
         quad_array.push_back(tmpquad);
}

	| Solution integer solname
{
#ifdef VERBOSE
	printf("YACC Solution number %d\n",$2);
#endif
	sol_head.sol_num = $2;

}
	| Component integer float  float  float float  integer float float
{
#ifdef VERBOSE
	printf("YACC Component \n");
#endif
	tmpcomponent.region = $2-1;
	tmpcomponent.xmin = $3;
	tmpcomponent.xmax = $4;
	tmpcomponent.ymax = -$5; /*invert y axis*/
	tmpcomponent.ymin = -$6; /*invert y axis*/
	tmpcomponent.direction  = $7;
	tmpcomponent.mole_begin = $8;
	tmpcomponent.mole_ratio = $9;
	component_array.push_back(tmpcomponent);
}

	| NodalSolution integer String value
{
#ifdef VERBOSE
	printf("YACC NodalSolution %d %s\n",$2,$3);
#endif
	tmpdata.index = $2-1;
	tmpdata.material = $3;
	sol_data.push_back(tmpdata);
	tmpdata.data_array.clear();
}
	| DataString String
{
#ifdef VERBOSE
	printf("YACC Transient Solution name %s\n",$2);
#endif
	strcpy(tran_sol.sol_name,$2);

}
	| Dataname namestring
{
#ifdef VERBOSE
	printf("YACC Transient Solution data name \n");
#endif

}
	| DataUnit unitstring
{
#ifdef VERBOSE
	printf("YACC Transient Solution data unit \n");
#endif
}
	| DataArray datalist
{
#ifdef VERBOSE
	printf("YACC Transient Solution data value \n");
#endif
}

	| MAT integer  parameters
{
#ifdef VERBOSE
	printf("YACC PARAMETER \n");
#endif
  	tmpparameter.region = $2;
  	parameter_array.push_back(tmpparameter);
  	tmpparameter.parameter_name_array.clear();
	tmpparameter.parameter_value_array.clear();
}
	| Ignore {}
	;


%%

int yyerror(char *s)
{
   printf("\nline %d unrecognized chars %s \n",yylineno,yylval.sval);
   return 0;
}

