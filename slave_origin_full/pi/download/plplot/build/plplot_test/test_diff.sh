#!/bin/bash
# Test suite to compare C examples with other language bindings
#
# Copyright (C) 2008 Andrew Ross
#
# This file is part of PLplot.
#
# PLplot is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published
# by the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# PLplot is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with PLplot; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

# Default device is psc. Will also work with svg.
device=${DEVICE:-psc}

usage()
{
echo '
Usage: test_diff.sh [OPTIONS]

Options:
   [--device=DEVICE] (DEVICE = any cmake-enabled device.  psc is the default)
   [--help]

Environment variables:
   DEVICE can be used to specify the device. This environment variable is
   overwritten by the --device option.
'
   exit $1
}

# Figure out what script options were specified by the user.

while test $# -gt 0; do
   if [ "ON" = "ON" ] ; then
      case "$1" in
      -*=*) optarg=${1#*=} ;;
      *) optarg= ;;
      esac
   else
      case "$1" in
      -*=*) optarg=`echo "$1" | sed 's/[-_a-zA-Z0-9]*=//'` ;;
      *) optarg= ;;
      esac
   fi
   case $1 in
      --device=*)
         device=$optarg
         ;;
      --help)
         usage 0 1>&2
         ;;
      *)
         usage 1 1>&2
         ;;
   esac
   shift
done


ret=0
# Comparison C results have no xsuffix.
xsuffix_c=

echo -e "\nComparison test using ${device} device\n"

if [ "${device}" = "svg" ] ; then
    usefam="yes"
    firstpage="01"
else
    usefam="no"
    firstpage=""
fi

# Compare C results with c++, fortran, java, octave, python, tcl, perl,
# ada, ocaml, lua and d results
for lang in c++ fortran java octave python tcl perl adastandard adatraditional ocaml lua d plrender; do
  # Check which suffix is used for this binding
  case $lang in
	c++)
	    xsuffix=
	    suffix=cxx
	    ;;
	fortran)
	    xsuffix=
	    suffix=f
	    ;;
	java)
	    xsuffix=
	    suffix=j
	    ;;
	octave)
	    xsuffix=
	    suffix=o
	    ;;
	python)
	    xsuffix=
	    suffix=p
	    ;;
	tcl)
	    xsuffix=
	    suffix=t
	    ;;
	perl)
	    xsuffix=
	    suffix=pdl
	    ;;
	adastandard)
	    xsuffix=standard
	    suffix=a
	    ;;
	adatraditional)
	    xsuffix=traditional
	    suffix=a
	    ;;
	ocaml)
	    xsuffix=
	    suffix=ocaml
	    ;;
	lua)
	    xsuffix=
	    suffix=lua
	    ;;
	d)
	    xsuffix=
	    suffix=d
	    ;;
	plrender)
	    xsuffix=
	    suffix=plm
	    ;;
    esac

    missing=""
    different=""
    diffstdout=""
    missingstdout=""

    # List of standard examples
    INDEX_LIST="00 01 02 03 04 05 06 07 08 09 10 11 12 13 14 14a 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 33"
    # Check if any examples exist for this language.
    EXAMPLES_EXIST="no"
    for index in ${INDEX_LIST} ; do
	  if [ -f x${xsuffix}${index}${suffix}${firstpage}.${device} ] ; then
	    EXAMPLES_EXIST="yes"
	    break
	  fi
    done
    if [ "$EXAMPLES_EXIST" = "yes" ] ; then
	  for index in ${INDEX_LIST} ; do
	    if [ ! -f x${xsuffix_c}${index}c${firstpage}.${device} ] ; then
                 echo "C example ${index} is missing"
	      else
		  if [ ! -f x${xsuffix}${index}${suffix}${firstpage}.${device} ] ; then
		    missing="${missing} ${index}"
		  else
                      if [ "${usefam}" = "yes" ] ; then
                          let i=1
                          p="01"
                          pages=""
                          while [ -f "x${xsuffix_c}${index}c${p}.${device}" ] ; do
                              pages="${pages} ${p}"
                              let i=i+1
                              printf -v p "%02d" $i
                          done
                      else
                          pages="xxx"
                      fi
                      isdiff="no"
                      for p in $pages ; do
                      if [ ${p} = "xxx" ] ; then
                          p=""
                      fi
		      if [ "ON" = "ON" ] ; then
                          if [ "${device}" = "psc" ] || [ "${device}" = "ps" ] ; then
		              # Skip first 190 bytes of comparison to ignore date stamp.
			      /usr/bin/cmp -s -i 190 x${xsuffix_c}${index}c${p}.${device} x${xsuffix}${index}${suffix}${p}.${device}
                          else
			      /usr/bin/cmp -s x${xsuffix_c}${index}c${p}.${device} x${xsuffix}${index}${suffix}${p}.${device}
                          fi
			  if [ $? != 0 ] && [ ! "${isdiff}" = "yes" ] ; then
			    different="${different} ${index}"
                            isdiff="yes"
			  fi
		    else
                          if [ "${device}" = "psc" ] || [ "${device}" = "ps" ] ; then
		              # Drop first 8 lines from comparison to ignore date stamp.
			      /usr/bin/tail -n +9 x${xsuffix_c}${index}c${p}.${device} > test1.${device}
			      /usr/bin/tail -n +9 x${xsuffix}${index}${suffix}${p}.${device} > test2.${device}
			      /usr/bin/diff -q test1.psc test2.psc 2>&1 > /dev/null
                          else
			      /usr/bin/diff -q x${xsuffix_c}${index}c${p}.${device} x${xsuffix}${index}${suffix}${p}.${device} 2>&1 > /dev/null
                          fi
			  if [ $? != 0 ] && [ ! "${isdiff}" = "yes" ] ; then
			    different="${different} ${index}"
                            isdiff="yes"
			  fi
		    fi
                    done
		    if [ "$index" != "14a" ] ; then
			  if [ -f x${xsuffix}${index}${suffix}_${device}.txt ] ; then
			    /usr/bin/diff -q x${xsuffix_c}${index}c_${device}.txt x${xsuffix}${index}${suffix}_${device}.txt 2>&1 > /dev/null
			    if [ $? != 0 ] ; then
				  diffstdout="${diffstdout} ${index}"
			    fi
			  else
			    missingstdout="${missingstdout} ${index}"
			  fi
		    fi
		  fi
	    fi
	  done
	echo "${lang}"
	echo "  Missing examples            : ${missing}"
	echo "  Differing graphical output  : ${different}"
	echo "  Missing stdout              : ${missingstdout}"
	echo "  Differing stdout            : ${diffstdout}"
	if [ "${different}" != "" -o "${diffstdout}" != "" ] ; then
	    ret=1
	fi
  fi
done

if [ "${ret}" != "0" ] ; then
    echo "WARNING: Some graphical or stdout results were different"
fi
