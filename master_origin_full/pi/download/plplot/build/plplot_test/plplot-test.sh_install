#!/bin/bash
# -*- mode: shell-script -*-
#
# Copyright (C) 2004-2017 Alan W. Irwin
# Copyright (C) 2004  Rafael Laboissiere
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

# test suite for executing all configured demos (either in plplot/tmp
# or else in the install area) capable of file output.  Our ultimate aim
# is to compare the generated files with results generated on a system
# that is known to give good results.  However, single-precision
# contamination (even when double is configured) and the consequent large
# roundoff error is currently hampering our cross-platform comparisons.
# So for now, this is only a test suite in the sense that it checks the
# files corresponding to our demo plots can all be generated.  It is
# up to the user to display those files with the appropriate software
# (e.g., gv for postscript files and kview or a browser for png, gif,
# or jpeg files on Linux systems) to make sure they display properly.

version=5.13.0

EXAMPLES_DIR=${EXAMPLES_DIR:-.}
SRC_EXAMPLES_DIR=${SRC_EXAMPLES_DIR:-.}
OUTPUT_DIR=${OUTPUT_DIR:-.}
device=${DEVICE:-psc}
# Check if we are in the build tree
if [ "`pwd`/.." -ef "/home/pi/download/plplot/build" ] ; then
  BUILD_TREE="YES"
else
  BUILD_TREE="NO"
fi
export EXAMPLES_DIR SRC_EXAMPLES_DIR OUTPUT_DIR device DEBUG_CMD BUILD_TREE

usage()
{
echo '
Usage: plplot-test.sh [OPTIONS]

Options:
   [--device=DEVICE] (DEVICE = any cmake-enabled device.  psc is the default)
   [--front-end=FE]  (FE = one of c, cxx, fortran java, octave, python, tcl, pdl, ada, ocaml, lua, d or plrender)
                     If this option is not specified, then all front-ends will
                     be tested.  More than one front-end may be given, like
                     this --front-end="c cxx"
   [--examples-dir=/path/to/examples/dir]
                     Specify path to the examples directory.  Defaults to "."
   [--src-examples-dir=/path/to/src/examples/dir]
                     Specify path to the src examples directory.  Defaults to "."
   [--output-dir=/path/to/output/dir]
                     Specify location where the resulting files are stored.
                     Defaults to "."
   [--interactive]   Run subset of C examples for interactive devices only.
   [--interactive_octave]   Run interactive octave examples for interactive devices only.
   [--verbose]	     Echo each PLplot example that is executed.
   [--debug="debug command"]
                     Run examples with given debug command.
   [--debug]         Run examples with default debug command (valgrind).
   [--version]
   [--help]

Environment variables:
   DEVICE, FRONT_END, EXAMPLES_DIR, SRC_EXAMPLES_DIR, and OUTPUT_DIR can be
   used to specify the devices, front-ends, the examples directory, the
   source examples directory (used for the special case of build-tree checks
   when that tree is separated from the source tree), and output directory.
   These environment variables are overridden by the options --device,
   --front-end, --examples-dir, --src-examples-dir, and --output-dir.
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
      --interactive)
         interactive=on
         ;;
      --interactive_octave)
         interactive_octave=on
         ;;
      --verbose)
         export verbose_test=on
         ;;
      --version)
         echo $version
	 exit 0
         ;;
      --front-end=*)
         FRONT_END=$optarg
         for i in $FRONT_END ; do
           [ $i = "c"         \
             -o $i = "cxx"    \
             -o $i = "fortran"    \
             -o $i = "java"   \
             -o $i = "octave" \
             -o $i = "python" \
             -o $i = "tcl"    \
             -o $i = "pdl"    \
             -o $i = "ada"    \
             -o $i = "ocaml"    \
             -o $i = "lua"    \
             -o $i = "d"    \
             -o $i = "plrender"    \
           ] \
           || usage 0 1>&2
         done
	 ;;
      --examples-dir=*)
         EXAMPLES_DIR=$optarg
         ;;
      --src-examples-dir=*)
         SRC_EXAMPLES_DIR=$optarg
         ;;
      --output-dir=*)
         OUTPUT_DIR=$optarg
         ;;
      --debug=*)
         DEBUG_CMD=$optarg
         ;;
      --debug)
         DEBUG_CMD="valgrind --leak-check=full --show-reachable=yes --log-file=valgrind.x%q{index}%q{lang}.%q{dsuffix}.log"
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

# This script is only designed to work when EXAMPLES_DIR is a directory
# with a subdirectory called "c".  Check whether this conditions is true.

if [ ! -d "$EXAMPLES_DIR"/c ] ; then
echo '
This script is only designed to work when the EXAMPLES_DIR environment
variable (overridden by option --examples-dir) is a directory with a
subdirectory called "c".  This condition has been violated.
'
exit 1
fi

# Find where the front-end scripts are by looking at the directory name of the
# current script.

if [ "" = "1" ] ; then
   scripts_dir=${0%/*}
else
   scripts_dir=`echo "$0" | sed 's:/[^/][^/]*$::'`
fi

if [ "$interactive" = "on" -o "$interactive_octave" = "on" ] ; then
    # List of interactive devices for PLplot that _might_ be enabled.
    PLD_aqt=OFF
    PLD_qtwidget=OFF
    PLD_ntk=OFF
    PLD_tk=OFF
    PLD_wingcc=OFF
    PLD_wincairo=OFF
    PLD_wxwidgets=OFF
    PLD_xcairo=OFF
    PLD_xwin=ON
	PLD_wingdi=OFF

    eval pld_device='$'PLD_$device
    if [ -z "$pld_device" ] ; then
	echo '
Never heard of an interactive device called '"$device"'.  Either this
is not a legitimate interactive device for PLplot or else
plplot-test.sh.cmake needs some maintenance to include this
interactive device in the list of possible PLplot interactive devices.
'
	exit 1
    fi

    if [ ! "$pld_device" = "ON" ] ; then
	echo '
PLD_'"$device"' is defined as '"$pld_device"'.  It must be ON (i.e., enabled
by your cmake configuration and built properly) before you can use this
script with DEVICE='"$device"'.
'
	exit 1
    fi

    status=0
    if [ "$interactive" = "on" ] ; then
	export cdir="$EXAMPLES_DIR"/c
	echo "Testing subset of C examples for device $device"
	script="$scripts_dir"/test_c_interactive.sh
    else
	# This case must be "$interactive_octave" = "on"
	export options=
	export octave=
        if [ "$BUILD_TREE" = "YES" ] ; then
	    export octavedir=\
"$SRC_EXAMPLES_DIR"/../bindings/octave/PLplot:\
"$SRC_EXAMPLES_DIR"/../bindings/octave/PLplot/support:\
"$SRC_EXAMPLES_DIR"/../bindings/octave/misc:\
"$SRC_EXAMPLES_DIR"/octave:\
"$EXAMPLES_DIR"/../bindings/octave:\
"$EXAMPLES_DIR"/../bindings/octave/PLplot
            if [ "OFF" = "ON" ] ; then
		octavedir=$octavedir:"$EXAMPLES_DIR"/../dll
	    fi
        else
            export octavedir=\
"":\
""/support:\
"":\
"":\
"$SRC_EXAMPLES_DIR"/octave
            if [ "OFF" = "ON" ] ; then
		octavedir=$octavedir:"/home/pi/plplot/install_directory/bin"
	    fi
        fi
	echo "Testing interactive octave examples for device $device"
	script="$scripts_dir"/test_octave_interactive.sh
    fi
    /bin/bash "$script" || status=1

    exit $status
fi

# These variables set by default assuming you are going to run this
# script from the installed demos directory $prefix/lib/plplot$version/examples.
cdir="$EXAMPLES_DIR"/c$VC_CTEST_DIRECTORY
cxxdir="$EXAMPLES_DIR"/c++$VC_CTEST_DIRECTORY
if [ "OFF" = "ON" ] ; then
    fortrandir="$EXAMPLES_DIR"/old_fortran
else
    fortrandir="$EXAMPLES_DIR"/fortran
fi
pythondir="$SRC_EXAMPLES_DIR"/python
tclexamplesdir="$EXAMPLES_DIR"/tcl
tcldir="$SRC_EXAMPLES_DIR"/tcl
perldir="$SRC_EXAMPLES_DIR"/perl
javadir="$EXAMPLES_DIR"/java
adadir="$EXAMPLES_DIR"/ada
ocamldir="$EXAMPLES_DIR"/ocaml
luadir="$SRC_EXAMPLES_DIR"/lua
ddir="$EXAMPLES_DIR"/d
octave=
if [ "$BUILD_TREE" = "YES" ] ; then
    octavedir=\
"$SRC_EXAMPLES_DIR"/../bindings/octave/PLplot:\
"$SRC_EXAMPLES_DIR"/../bindings/octave/PLplot/support:\
"$SRC_EXAMPLES_DIR"/../bindings/octave/misc:\
"$SRC_EXAMPLES_DIR"/octave:\
"$EXAMPLES_DIR"/../bindings/octave:\
"$EXAMPLES_DIR"/../bindings/octave/PLplot
    if [ "OFF" = "ON" ] ; then
	octavedir=$octavedir:"$EXAMPLES_DIR"/../dll
    fi
else
    export octavedir=\
"":\
""/support:\
"":\
"":\
"$SRC_EXAMPLES_DIR"/octave
    if [ "OFF" = "ON" ] ; then
	octavedir=$octavedir:"/home/pi/plplot/install_directory/bin"
    fi
fi
export cdir cxxdir fortrandir pythondir javadir octave octavedir tcldir perldir adadir ocamldir luadir ddir

fe=""

# List of non-interactive (i.e., file) devices for PLplot that
# _might_ be enabled.  For completeness you may want to specify all devices
# here, but be sure to comment out the interactive ones since they are
# handled by the above --interactive logic.

#interactive PLD_aqt=OFF
PLD_cgm=OFF
PLD_epsqt=OFF
PLD_pdfqt=OFF
#interactive PLD_qtwidget=OFF
PLD_bmpqt=OFF
PLD_jpgqt=OFF
PLD_pngqt=OFF
PLD_ppmqt=OFF
PLD_tiffqt=OFF
PLD_svgqt=OFF
PLD_gif=OFF
PLD_jpeg=OFF
#not a file device PLD_mem=ON
#not a file device PLD_memcairo=OFF
#interactive PLD_ntk=OFF
PLD_null=ON
PLD_pdf=OFF
PLD_pdfcairo=OFF
PLD_plmeta=OFF
PLD_png=OFF
PLD_pngcairo=OFF
PLD_ps=ON
# special case
PLD_psc=ON
PLD_pscairo=OFF
PLD_epscairo=OFF
PLD_pstex=OFF
PLD_psttf=OFF
# special case
PLD_psttfc=OFF
PLD_svg=ON
PLD_svgcairo=OFF
#interactive PLD_tk=OFF
#interactive PLD_tkwin=OFF
#interactive PLD_wincairo=OFF
#interactive PLD_wingcc=OFF
PLD_wxpng=OFF
#interactive PLD_wxwidgets=OFF
#interactive PLD_xcairo=OFF
PLD_xfig=ON
#interactive PLD_xwin=ON

eval pld_device='$'PLD_$device
if [ -z "$pld_device" ] ; then
echo '
Never heard of a file device called '"$device"'.  Either this is not a
legitimate file (i.e. non-interactive) device for PLplot or else
plplot-test.sh.cmake needs some maintenance to include this file device in
the list of possible PLplot file devices.
'
exit 1
fi

if [ ! "$pld_device" = "ON" ] ; then
echo '
PLD_'"$device"' is defined as '"$pld_device"'.  It must be ON (i.e., enabled
by your cmake configuration and built properly) before you can use this
script with DEVICE='"$device"'.
'
exit 1
fi

# Some devices require special options others do not.
case "$device" in
   png|pngcairo|epscairo|jpeg|xfig|svg|svgcairo|bmpqt|jpgqt|pngqt|ppmqt|tiffqt|svgqt|epsqt|pdfqt)
      options="-fam -fflen 2"
      ;;
   gif)
   # gif standard is limited to 256 colours so cannot have aliasing turned
   # on (which uses additional colours) for examples which already have
   # a substantial number of colours (such as example 2 with its 116 colours).
      options="-fam -fflen 2 -drvopt smooth=0"
      ;;
   *)
      options=
      ;;
esac
dsuffix=$device
export dsuffix options

# Find out what front-ends have been configured
if [ -z "$FRONT_END" ] ; then
   FRONT_END=c
   test "ON" = "ON" && FRONT_END="$FRONT_END cxx"
   test "ON" = "ON"    && FRONT_END="$FRONT_END fortran"
   test "OFF" = "ON"   && FRONT_END="$FRONT_END java"
   test "OFF" = "ON" && FRONT_END="$FRONT_END octave"
   test "OFF" = "ON" && FRONT_END="$FRONT_END python"
   test "OFF" = "ON"    && FRONT_END="$FRONT_END tcl"
   test "OFF" = "ON"    && FRONT_END="$FRONT_END pdl"
   test "OFF" = "ON"    && FRONT_END="$FRONT_END ada"
   test "OFF" = "ON"  && FRONT_END="$FRONT_END ocaml"
   test "OFF" = "ON"  && FRONT_END="$FRONT_END lua"
   test "OFF" = "ON"  && FRONT_END="$FRONT_END d"
   test "OFF" = "ON" && FRONT_END="$FRONT_END plrender"
fi

# Call the front-end scripts
status=0

if [ "" = "1" ] ; then
 exesuffix=".exe"
else
 exesuffix=
fi
critical_examples="14 17 29"
export critical_examples exesuffix
for i in $FRONT_END ; do
   echo "Testing front-end $i"
   script="$scripts_dir"/test_$i.sh
   /bin/bash "$script" || status=1
done

exit $status
