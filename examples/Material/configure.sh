#!/bin/bash
CONFIG_H="config.h"
MAKE_DEFS="make.defs"

cat>${CONFIG_H} <<_ACEOF
#ifndef CONFIG_H
#define CONFIG_H

_ACEOF

echo "" >${MAKE_DEFS}

# test processor
OSSTRING=`uname -s`
case "$OSSTRING" in
	i386|x86_64)
	echo "#define I386" >>${CONFIG_H};
	;;
	powerpc)
	echo "#define POWERPC" >>${CONFIG_H};
	;;
esac

# test OS
OSSTRING=`uname -s`
case "$OSSTRING" in
	Darwin)
	echo "#define DARWIN" >>${CONFIG_H};
	echo "object_libext = .o" >> ${MAKE_DEFS};
	echo "obj-suffix = .o" >> ${MAKE_DEFS};
	echo "static_libext = .a" >> ${MAKE_DEFS};
	echo "shared_libext = .so" >> ${MAKE_DEFS};
	echo "DLLFLAG = -fPIC" >> ${MAKE_DEFS} ;
	echo "CXXSHAREDFLAG  = -dynamiclib -undefined dynamic_lookup -single_module" >> ${MAKE_DEFS} ;
	echo "EXE = " >> ${MAKE_DEFS};
	;;
	Linux)
	echo "#define LINUX" >>${CONFIG_H};
	echo "objext = .o" >> ${MAKE_DEFS};
	echo "obj-suffix = .o" >> ${MAKE_DEFS};
	echo "static_libext = .a" >> ${MAKE_DEFS};
	echo "shared_libext = .so" >> ${MAKE_DEFS};
	echo "DLLFLAG = -fPIC" >> ${MAKE_DEFS} ;
	echo "CXXSHAREDFLAG  = -shared" >> ${MAKE_DEFS} ;
	echo "EXE = " >> ${MAKE_DEFS};
	;;
	CYGWIN*)
	echo "#define CYGWIN" >>${CONFIG_H};
	echo "object_libext = .o" >> ${MAKE_DEFS};
	echo "obj-suffix = .o" >> ${MAKE_DEFS};
	echo "static_libext = .lib" >> ${MAKE_DEFS};
	echo "shared_libext = .dll" >> ${MAKE_DEFS};
	;;
esac

# test Compiler
if [ "x$CC" == "x" ]
then
	echo "CC=gcc" >> ${MAKE_DEFS} ;
	echo "CXX=g++" >> ${MAKE_DEFS} ;
fi
echo "AR=ar" >> ${MAKE_DEFS} ;
echo "RM=rm" >> ${MAKE_DEFS} ;
echo "LINK=g++" >> ${MAKE_DEFS} ;

if [ "$METHOD" == "dbg" ]
then
	echo "CFLAGS= -g" >> ${MAKE_DEFS} ;
	echo "CXXFLAGS = -g" >> ${MAKE_DEFS} ;
else
	echo "CFLAGS= " >> ${MAKE_DEFS} ;
	echo "CXXFLAGS= " >> ${MAKE_DEFS} ;
fi



cat >>$CONFIG_H <<_ACCOF

#endif
_ACCOF


