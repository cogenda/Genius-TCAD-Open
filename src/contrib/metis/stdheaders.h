/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * stdheaders.h
 *
 * This file includes all necessary header files
 *
 * Started 8/27/94
 * George
 *
 * $Id: stdheaders.h,v 1.1 2008/04/01 10:22:19 gdiso Exp $
 */


#include <stdio.h>
#ifdef __STDC__
#include <stdlib.h>
#elif defined(__APPLE__)
#include <stdlib.h>
#else
#include <malloc.h>
#endif
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <stdarg.h>
#include <time.h>

