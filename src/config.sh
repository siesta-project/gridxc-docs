#!/bin/sh
#
##set -x
#
# Get absolute path of this script, as that will be the Src directory to use
# as reference when copying files.
# 
#
srcdir=$(
cd -P -- "$(dirname -- "$0")" &&
pwd -P
)
# The above construct is more robust than:  srcdir=$(dirname $0)
# (It will work if $0 is "../Src", since we want an *absolute* path
#
# Get also the absolute path of the object directory
#
objdir=$(
cd -P -- "$(pwd)" &&
pwd -P
)
#
destdir=$objdir
#
# Copy test material
#
(cd $srcdir; cd .. ; cp -rp Testers ${destdir} )
#
# Copy and change name of makefiles. The '.alt' suffix was used to avoid
# clashes with the auto-tools automatically-generated files.
#
(
 cp $srcdir/makefile.alt $destdir/makefile
 cd $destdir/Testers; mv makefile.alt makefile
)
#
# Replicate any .inc files
# This is needed in some systems with broken include file import heuristics
#
(cd $srcdir;
  for i in $(find . -name '*.inc' ); do
    relpath=${i%/*}
    mkdir -p ${destdir}/$relpath
    cp -f $relpath/*.inc ${destdir}/$relpath
  done
)
#
#
# Copy other needed top-level files
#
cp -p $srcdir/build.sh ${destdir}
cp -p $srcdir/gridxc.mk.in ${destdir}
cp -p $srcdir/libxc.mk ${destdir}
cp -p $srcdir/top.gridxc.mk.in ${destdir}
#
# Set the appropiate variables in the build makefile
#
sed "s#VPATH=\.#VPATH=${srcdir}#g" ${srcdir}/makefile.alt | \
sed "s#MAIN_OBJDIR=\.#MAIN_OBJDIR=${objdir}#g" > ${destdir}/makefile

#
echo " *** Compilation setup done. "

## for i in $(find . \( -path Tutorial -o -path Examples \) -prune -o \
