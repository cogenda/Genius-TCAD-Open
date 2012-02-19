#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>


typedef double (* FUNC)(double);

#define XPRECISION 1.0e-15
#define FPRECISION 1.0e-100
#define MAX_ITERATIONS 100

#ifndef MAXDOUBLE
#define MAXDOUBLE   1.797693E+308
#endif


char *apszHeader1[] = {
"/*--------------------------------------------------------------------------",
"",
" PROJECT:  SIMULATION GENERATION FRAMEWORK (SIMGEN)",
"",
" SUBSYSTEM:  brkpnts.exe",
" FILE:	      brkpnts.c",
" AUTHOR:     Kevin M. Kramer",
"",
" DESCRIPTION:",
"",
" This module contains the range limit constant for the Bernoulli, Aux1 and",
" Aux2 functions as well as their derivatives.	For more information consult",
" \"Analysis and Simulation of Electronic Devices\" by Siegfried Selberherr",
" pages 158 and 159.",
"",
"--------------------------------------------------------------------------*/",
"",
"",
"#ifndef __brkpnts_h",
"#define __brkpnts_h",
"",
NULL };


int sign(double x)
{
  if	  (x < 0.0) return(-1);
  else if (x > 0.0) return(+1);
  else		    return(0);
}


double Bisection(FUNC func1, FUNC func2, double Xpos, double Xneg)
{
  double Fpos = func1(Xpos) - func2(Xpos);
  double Fneg = func1(Xneg) - func2(Xneg);
  double Xmid, Fmid, Xlast;

  if	  (Fpos == 0.0) return(Xpos);
  else if (Fneg == 0.0) return(Xneg);
  else if ((Fpos > 0.0) && (Fneg < 0.0)) ;
  else if ((Fpos < 0.0) && (Fneg > 0.0))
  {
    Xmid = Xpos;
    Xpos = Xneg;
    Xneg = Xmid;
  }
  else
  {
    fprintf(stderr, "Initial interval may not contain a root\n");
    exit(1);
  }

  Xlast = 0.0;
  do
  {
    Xmid = 0.5 * (Xpos + Xneg);
    Fmid = func1(Xmid) - func2(Xmid);
    if	    (Fmid < 0.0)  Xneg = Xmid;
    else if (Fmid > 0.0)  Xpos = Xmid;
    if (Xlast == Xmid) return(Xmid);
    else Xlast = Xmid;
  } while (Xneg != Xpos);

  return(Xmid);

}


double Secant(FUNC func1, FUNC func2, double x1)
{
  double slope, dx, x3, f3;
  int	 s3, iteration;

  double x2 = 0.9 * x1;
  double f1 = func1(x1) - func2(x1);
  double f2 = func1(x2) - func2(x2);
  int	 s2 = sign(x2);

  for(;;)
  {
    iteration = 0;
    slope = (f2 - f1) / (x2 - x1);
    dx = f2 / slope;
    x3 = x2 - dx;
    f3 = func1(x3) - func2(x3);
    s3 = sign(x3);

    while ((fabs(f3) >= fabs(f2)) || (s3 != s2))
    {
      dx /= 1.97;
      x3 += dx;
      f3  = func1(x3) - func2(x3);
      s3  = sign(x3);

      if (++iteration > MAX_ITERATIONS)
      {
  	if (fabs(f2) <= 100.0 * XPRECISION)
	{
		return(x2);
	}
  	fprintf(stderr, "ERROR:  secant method not converging!\n");
  	exit(1);
      }

    }

    x1 = x2;
    x2 = x3;
    f1 = f2;
    f2 = f3;

    if ((fabs(dx / x2) <= XPRECISION) || (fabs(f2) <= XPRECISION)) break;
  }

  return(x3);

}


double Asymptotic(FUNC func1, FUNC func2, double x, double dx)
{
  while (1)
  {
    if (fabs(dx/x) <= XPRECISION) return(x);
    while ( fabs(func1(x) - func2(x)) >= FPRECISION)
      x += dx;

    dx *= -0.1;
    if (fabs(dx/x) <= XPRECISION) return(x);
    while (fabs(func1(x) - func2(x)) <= FPRECISION)
      x += dx;
    dx *= -0.1;
  }
}


double Bbp0a(double x)	{ return exp(x) - 1.0; }
double Bbp0b(double x)	{ return - 1.0; }

double Bbp1a(double x)	{ return x / (exp(x) - 1.0); }
double Bbp1b(double x)	{ return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0)); }

double Bbp2a(double x)	{ return 1.0 - x/2.0 * (1.0 - x/6.0 * (1.0 - x*x/60.0)); }
double Bbp2b(double x)	{ return x * exp(-x) / (1.0 - exp(-x)); }

double Bbp3a(double x)	{ return 1.0 - exp(-x); }
double Bbp3b(double x)	{ return 1.0; }

double Bbp4a(double x)	{ return x * exp(-x); }
double Bbp4b(double x)	{ return 0.0; }

double dBbp0a(double x) { return (1.0 - x) * exp(x) - 1.0; }
double dBbp0b(double x) { return -1.0; }

double dBbp2a(double x) { return ((1.0 - x) * exp(x) - 1.0) / ((exp(x) - 1.0) * (exp(x) - 1.0)); }
double dBbp2b(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }

double dBbp3a(double x) { return -0.5 + x/6.0 * (1.0 - x*x/30.0); }
double dBbp3b(double x) { return (exp(-x) * (1.0 - x) - exp(-2.0 * x)) / ((1.0 - exp(-x)) * (1.0 - exp(-x))); }

double dBbp5a(double x) { return exp(-x) * (1.0 - x) - exp(-2.0 * x); }
double dBbp5b(double x) { return 0.0; }

double AUX1bp0a(double x) { return x / sinh(x); }
double AUX1bp0b(double x) { return 1.0 - x*x/6.0 * (1.0 - 7.0*x*x/60.0); }

double dAUX1bp0a(double x) { return (sinh(x) - x*cosh(x)) / (sinh(x) * sinh(x)); }
double dAUX1bp0b(double x) { return -x/3.0 * (1.0 - 7*x*x/30.0); }

double dAUX1bp2a(double x) { return (sinh(x) - x*cosh(x)) / (sinh(x) * sinh(x)); }
double dAUX1bp2b(double x) { return -2.0 * (1 + x) * exp(x); }

double dAUX1bp3a(double x) { return (sinh(x) - x*cosh(x)) / (sinh(x) * sinh(x)); }
double dAUX1bp3b(double x) { return  2.0 * (1 - x) * exp(-x); }

double dAUX1bp4a(double x) { return -2.0 * (1 + x) * exp(x); }
double dAUX1bp4b(double x) { return 0.0; }

double dAUX1bp5a(double x) { return  2.0 * (1 - x) * exp(-x); }
double dAUX1bp5b(double x) { return 0.0; }

double AUX2bp0a(double x) { return 1.0; }
double AUX2bp0b(double x) { return 1.0 + exp(x); }

double AUX2bp1a(double x) { return 1.0 + exp(x); }
double AUX2bp1b(double x) { return exp(x); }

double AUX2bp2a(double x) { return exp(-x); }
double AUX2bp2b(double x) { return 0.0; }

double dAUX2bp0a(double x) { return exp(x); }
double dAUX2bp0b(double x) { return 0.0; }


int main ( int cArg, char *apszArg[] )
{
  FILE	 *file;
  char	 **ppsz;
  double B[5], dB[6], AUX1[2], dAUX1[6], AUX2[3], dAUX2[4];

  fprintf(stderr, "Bernoulli and Aux function limit evaluation, (c) 1992 K. Kramer\n");


  fprintf(stderr, "BP0_BERN: Asymptotic\n");
  B[0]	= Asymptotic(Bbp0a, Bbp0b,  0.0e+00, -1.0e+02);

  fprintf(stderr, "BP1_BERN: Secant\n");
  B[1]	= Secant    (Bbp1a, Bbp1b, -1.0e+00);

  fprintf(stderr, "BP2_BERN: Secant\n");
  B[2]	= Secant    (Bbp2a, Bbp2b, +1.0e+00);

  fprintf(stderr, "BP3_BERN: Asymptotic\n");
  B[3]	= Asymptotic(Bbp3a, Bbp3b,  0.0e+00, +1.0e+02);

  fprintf(stderr, "BP4_BERN: Asymptotic\n");
  B[4]	= Asymptotic(Bbp4a, Bbp4b,  1.0e+00, +1.0e+03);



  fprintf(stderr, "BP0_DBERN: Asymptotic\n");
  dB[0] = Asymptotic(dBbp0a, dBbp0b,  0.0e+00, -1.0e+02);

  dB[1] = B[0];

  fprintf(stderr, "BP2_DBERN: Secant\n");
  dB[2] = Secant    (dBbp2a, dBbp2b, -1.0e+00);
/*dB[2] = Bisection (dBbp2a, dBbp2b, dB[1], -1.0e-6); */

  fprintf(stderr, "BP3_DBERN: Secant\n");
  dB[3] = Secant    (dBbp3a, dBbp3b, +1.0e+00);

  dB[4] = B[3];

  fprintf(stderr, "BP5_DBERN: Asymptotic\n");
  dB[5] = Asymptotic(dBbp5a, dBbp5b,  1.0e+00, +1.0e+03);


  fprintf(stderr, "BP0_AUX1: Secant\n");
  AUX1[0] = Secant(AUX1bp0a, AUX1bp0b, -1.0e+00);

  fprintf(stderr, "BP1_AUX1: Secant\n");
  AUX1[1] = Secant(AUX1bp0a, AUX1bp0b,	1.0e+00);


  fprintf(stderr, "BP0_DAUX1: Secant\n");
  dAUX1[0] = Secant(dAUX1bp0a, dAUX1bp0b, -1.0e+00);

  fprintf(stderr, "BP1_DAUX1: Secant\n");
  dAUX1[1] = Secant(dAUX1bp0a, dAUX1bp0b,  1.0e+00);

  fprintf(stderr, "BP2_DAUX1: Asymptotic\n");
  dAUX1[2] = Asymptotic(dAUX1bp2a, dAUX1bp2b,  -1.0e+00, -1.0e+02);

  fprintf(stderr, "BP3_DAUX1: Asymptotic\n");
  dAUX1[3] = Asymptotic(dAUX1bp3a, dAUX1bp3b,  1.0e+00,  1.0e+02);

  fprintf(stderr, "BP4_DAUX1: Asymptotic\n");
  dAUX1[4] = Asymptotic(dAUX1bp4a, dAUX1bp4b,  0.0e+00, -1.0e+02);

  fprintf(stderr, "BP5_DAUX1: Asymptotic\n");
  dAUX1[5] = Asymptotic(dAUX1bp5a, dAUX1bp5b,  0.0e+00,  1.0e+02);

  fprintf(stderr, "BP0_AUX2: Asymptotic\n");
  AUX2[0] = Asymptotic(AUX2bp0a, AUX2bp0b,  0.0e+00, -1.0e+02);

  fprintf(stderr, "BP1_AUX2: Asymptotic\n");
  AUX2[1] = Asymptotic(AUX2bp1a, AUX2bp1b,  0.0e+00,  1.0e+02);

  fprintf(stderr, "BP2_AUX2: Asymptotic\n");
  AUX2[2] = Asymptotic(AUX2bp2a, AUX2bp2b,  0.0e+00,  1.0e+02);


  fprintf(stderr, "BP0_DAUX2: Asymptotic\n");
  dAUX2[0] = Asymptotic(dAUX2bp0a, dAUX2bp0b,  0.0e+00, -1.0e+02);

  dAUX2[1] = AUX2[0];

  dAUX2[2] = AUX2[1];

  dAUX2[3] = AUX2[2];


  fprintf(stderr, "Test completed, writing to output file.\n");


  /* write header file */
  file = fopen("brkpnts.h", "wt");
  for(ppsz = apszHeader1; *ppsz; ++ppsz) fprintf(file, "%s\n", *ppsz);
  fprintf(file, "#define BP0_BERN    %.15le\n", B[0]);
  fprintf(file, "#define BP1_BERN    %.15le\n", B[1]);
  fprintf(file, "#define BP2_BERN    %.15le\n", B[2]);
  fprintf(file, "#define BP3_BERN    %.15le\n", B[3]);
  fprintf(file, "#define BP4_BERN    %.15le\n", B[4]);
  fprintf(file, "#define BP0_DBERN   %.15le\n", dB[0]);
  fprintf(file, "#define BP1_DBERN   %.15le\n", dB[1]);
  fprintf(file, "#define BP2_DBERN   %.15le\n", dB[2]);
  fprintf(file, "#define BP3_DBERN   %.15le\n", dB[3]);
  fprintf(file, "#define BP4_DBERN   %.15le\n", dB[4]);
  fprintf(file, "#define BP5_DBERN   %.15le\n", dB[5]);
  fprintf(file, "#define BP0_AUX1    %.15le\n", AUX1[0]);
  fprintf(file, "#define BP1_AUX1    %.15le\n", AUX1[1]);
  fprintf(file, "#define BP0_DAUX1   %.15le\n", dAUX1[0]);
  fprintf(file, "#define BP1_DAUX1   %.15le\n", dAUX1[1]);
  fprintf(file, "#define BP2_DAUX1   %.15le\n", dAUX1[2]);
  fprintf(file, "#define BP3_DAUX1   %.15le\n", dAUX1[3]);
  fprintf(file, "#define BP4_DAUX1   %.15le\n", dAUX1[4]);
  fprintf(file, "#define BP5_DAUX1   %.15le\n", dAUX1[5]);
  fprintf(file, "#define BP0_AUX2    %.15le\n", AUX2[0]);
  fprintf(file, "#define BP1_AUX2    %.15le\n", AUX2[1]);
  fprintf(file, "#define BP2_AUX2    %.15le\n", AUX2[2]);
  fprintf(file, "#define BP0_DAUX2   %.15le\n", dAUX2[0]);
  fprintf(file, "#define BP1_DAUX2   %.15le\n", dAUX2[1]);
  fprintf(file, "#define BP2_DAUX2   %.15le\n", dAUX2[2]);
  fprintf(file, "#define BP3_DAUX2   %.15le\n", dAUX2[3]);
  fprintf(file, "#define BP0_MISC    %.15le\n", log(MAXDOUBLE));
  fprintf(file, "\n#endif\n");
  fclose(file);


  /* write text file */
  file = fopen("brkpnts.txt", "wt");
  fprintf(file, "This program determines how to evaluate the following functions to machine\n");
  fprintf(file, "precision.  The breakpoints of these functions are dependent upon the host\n");
  fprintf(file, "computer's floating point architecture and floating point library.\n");
  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "		       x\n");
  fprintf(file, "	     B(x) = -------\n");
  fprintf(file, "		    e^x - 1\n");
  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "	     d	      (1-x)*e^x - 1\n");
  fprintf(file, "	     --B(x) = -------------\n");
  fprintf(file, "	     dx	       (e^x - 1)^2\n");
  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "			   x\n");
  fprintf(file, "	     Aux1(x) =	-------\n");
  fprintf(file, "			sinh(x)\n");
  fprintf(file, "\n");
  fprintf(file, "	     d		 sinh(x) - x*cosh(x)\n");
  fprintf(file, "	     --Aux1(x) = -------------------\n");
  fprintf(file, "	     dx		     (sinh(x))^2\n");
  fprintf(file, "\n");
  fprintf(file, "			  1\n");
  fprintf(file, "	     Aux2(x) = -------\n");
  fprintf(file, "		       1 + e^x\n");
  fprintf(file, "\n");
  fprintf(file, "	     d		   - e^x\n");
  fprintf(file, "	     --Aux2(x) = -----------\n");
  fprintf(file, "	     dx		 (1 + e^x)^2\n");
  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "To achieve machine precision, the above functions should ");
  fprintf(file, " be evaluated as\nfollows:\n");
  fprintf(file, "\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "	      |	 -x						 x <= %+9.1le\n", B[0]);
  fprintf(file, "	      |	 x / (e^x - 1)			    %+9.1le <  x <  %+9.1le\n", B[0], B[1]);
  fprintf(file, "     B(x) = <	 1 - x/2*(1 - x/6*(1 - x*x/60))	    %+9.1le <= x <= %+9.1le\n", B[1], B[2]);
  fprintf(file, "	      |	 x*e^-x / (1 - e^-x)		    %+9.1le <  x <  %+9.1le\n", B[2], B[3]);
  fprintf(file, "	      |	 x*e^-x				    %+9.1le <= x <  %+9.1le\n", B[3], B[4]);
  fprintf(file, "	      |	 0				    %+9.1le <= x\n", B[4]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "	      |	 -1						 x <= %+9.1le\n", dB[0]);
  fprintf(file, "	      |	 {(1-x)*e^x - 1}		    %+9.1le <  x <= %+9.1le\n", dB[0], dB[1]);
  fprintf(file, "   d	      |	 {(1-x)*e^x - 1} / (e^x - 1)^2	    %+9.1le <  x <  %+9.1le\n", dB[1], dB[2]);
  fprintf(file, "   --B(x) = <	 -1/2 + x/6*(1 - x*x/30)	    %+9.1le <= x <= %+9.1le\n", dB[2], dB[3]);
  fprintf(file, "   dx	      |	 {(1-x)*e^-x - e^-2x}/(1 - e^-x)^2  %+9.1le <  x <  %+9.1le\n", dB[3], dB[4]);
  fprintf(file, "	      |	 {(1-x)*e^-x - e^-2x}		    %+9.1le <= x <  %+9.1le\n", dB[4], dB[5]);
  fprintf(file, "	      |	 0				    %+9.1le <= x\n", dB[5]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "	      |	 x / sinh(x)					 x <= %+9.1le\n", AUX1[0]);
  fprintf(file, "  Aux1(x) = <	 1 - x*x/6*(1 - 7*x*x/60)	    %+9.1le <  x <  %+9.1le\n", AUX1[0], AUX1[1]);
  fprintf(file, "	      |	 x / sinh(x)			    %+9.1le <= x\n", AUX1[1]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "d	      |	 {sinh(x) - x*cosh(x)}/{sinh(x)}^2		 x <= %+9.1le\n", dAUX1[0]);
  fprintf(file, "--Aux1(x) = <	 -x/3*(1 - 7*x*x/30)		    %+9.1le <  x <  %+9.1le\n", dAUX1[0], dAUX1[1]);
  fprintf(file, "dx	      |	 {sinh(x) - x*cosh(x)}/{sinh(x)}^2  %+9.1le <= x\n", dAUX1[1]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "	      |	 1						 x <= %+9.1le\n", AUX2[0]);
  fprintf(file, "  Aux2(x) = <	 1 / (1 + e^x)			    %+9.1le <  x <  %+9.1le\n", AUX2[0], AUX2[1]);
  fprintf(file, "	      |	 e^-x				    %+9.1le <  x <  %+9.1le\n", AUX2[1], AUX2[2]);
  fprintf(file, "	      |	 0				    %+9.1le <= x\n", AUX2[2]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fprintf(file, "	       /\n");
  fprintf(file, "	      |	 0						 x <= %+9.1le\n", dAUX2[0]);
  fprintf(file, "d	      |	 - e^x				    %+9.1le <  x <  %+9.1le\n", dAUX2[0], dAUX2[1]);
  fprintf(file, "--Aux2(x) = <	 - e^x / {1 + e^x}^2		    %+9.1le <= x <= %+9.1le\n", dAUX2[1], dAUX2[2]);
  fprintf(file, "dx	      |	 - e^-x				    %+9.1le <  x <  %+9.1le\n", dAUX2[2], dAUX2[3]);
  fprintf(file, "	      |	 0				    %+9.1le <= x\n", dAUX2[3]);
  fprintf(file, "	       \\\n");
  fprintf(file, "\n");
  fclose(file);

  return(0);
}

