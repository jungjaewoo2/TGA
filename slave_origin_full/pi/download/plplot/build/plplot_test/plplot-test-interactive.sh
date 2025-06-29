#!/bin/bash
# -*- mode: shell-script -*-
# Copyright (C) 2009-2017 Alan W. Irwin
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

# Test suite of PLplot interactive stuff that cannot be tested with
# file output device drivers.

EXAMPLES_DIR=${EXAMPLES_DIR:-.}
SRC_EXAMPLES_DIR=${SRC_EXAMPLES_DIR:-.}

# OPTIONAL tests.
TEST_PLTCL_STANDARD_EXAMPLES=1
TEST_TCLSH_STANDARD_EXAMPLES=1
TEST_PLSERVER_STANDARD_EXAMPLES=1
TEST_WISH_STANDARD_EXAMPLES=0
TEST_PLSERVER_RUNALLDEMOS=0
TEST_WISH_RUNALLDEMOS=0
TEST_QT_EXAMPLE=0

usage()
{
echo '
Usage: plplot-test-interactive.sh [OPTIONS]

Options:
   [--device=DEVICES_LIST]
       where DEVICES_LIST is one of more devices specified as a
       blank-delimited string (e.g., --device=xwin or
       --device="qtwidget extqt").  Note each specified device must be
       taken from the following list of eligible interactive devices:

       xwin, tk, xcairo, wxwidgets, qtwidget, extcairo or extqt.

   [--help]

Environment variables:
   DEVICES_LIST can be used to specify the device(s).
   This environment variable is overridden by the option --device.

N.B. All members of DEVICES_LIST (whether specified by the
DEVICES_LIST environment variable or by the --device option _must_ be
configured.  If neither the environment variable or --device option is
specified, then every _configured_ device from the above list is used.
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
         DEVICES_LIST=$optarg
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

# Set up interactive and/or external devices that are used for tests.
# Cannot use loop for this because of configuration.

PLD_xwin=ON
PLD_tk=OFF
PLD_ntk=OFF
PLD_xcairo=OFF
PLD_wxwidgets=OFF
PLD_qtwidget=OFF
PLD_extcairo=OFF
PLD_extqt=OFF

# MAINTENANCE: These blank-delimited strings must be consistent with
# previous configured list of devices.
POSSIBLE_INTERACTIVE_DEVICES_LIST="xwin tk ntk xcairo qtwidget"
if [ "1" -ne 0 ] ; then
    POSSIBLE_INTERACTIVE_DEVICES_LIST="$POSSIBLE_INTERACTIVE_DEVICES_LIST wxwidgets"
fi
POSSIBLE_DEVICES_LIST="$POSSIBLE_INTERACTIVE_DEVICES_LIST extcairo extqt"

# Default DEVICES_LIST is all eligible devices if environment variable
# not specified and --devices option not specified.
if [ -z "$DEVICES_LIST" ] ; then
    DEVICES_LIST=
    for device in $POSSIBLE_DEVICES_LIST ; do
	eval pld_device='$'PLD_$device
	test "$pld_device" = "ON" && DEVICES_LIST="$DEVICES_LIST $device"
    done
fi

# Check that everything in DEVICES_LIST is a configured  device.
for device in $DEVICES_LIST ; do
    eval pld_device='$'PLD_$device
    if [ ! "$pld_device" = "ON" ] ; then
	echo "$device is either not configured or not used for interactive tests."
	usage 1 1>&2
    fi
done

# Turn off all devices not mentioned in DEVICES_LIST.
for device in $POSSIBLE_DEVICES_LIST ; do
    eval "PLD_$device=OFF"
done

for device in $DEVICES_LIST ; do
    eval "PLD_$device=ON"
done

INTERACTIVE_DEVICES_LIST=
for device in $POSSIBLE_INTERACTIVE_DEVICES_LIST ; do
    eval pld_device='$'PLD_$device
    test "$pld_device" = "ON" && \
	INTERACTIVE_DEVICES_LIST="$INTERACTIVE_DEVICES_LIST $device"
done

OVERALL_STATUS_CODE=0
for device in $INTERACTIVE_DEVICES_LIST ; do
    ./plplot-test.sh --verbose --interactive --device=$device
    OVERALL_STATUS_CODE=$?
done

INDEX_LIST=
COUNT=0
if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/c"
    COMMAND[$COUNT]="./extXdrawable_demo"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/c"
    COMMAND[$COUNT]="./ext-cairo-test"
    COUNT=$(( $COUNT + 1 ))
fi

# ENABLE_wxwidgets is ON only if ENABLE_cxx is ON.
if [ "OFF" = "ON" -a "1" -ne 0 ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/c++"
    COMMAND[$COUNT]="./"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" -a "ON" = "ON" -a "#" != "#" -a "$TEST_QT_EXAMPLE" -ne 0 ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/c++"
    COMMAND[$COUNT]="./qt_example"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/tk"
    COMMAND[$COUNT]="./xtk01 -f tk01"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/tk"
    COMMAND[$COUNT]="./xtk02 -f tk02"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/tk"
    COMMAND[$COUNT]="plserver -f tk03"
    COUNT=$(( $COUNT + 1 ))
fi

if [ "OFF" = "ON" ] ; then
    INDEX_LIST="$INDEX_LIST $COUNT"
    DIRECTORY[$COUNT]="${EXAMPLES_DIR}/tk"
    COMMAND[$COUNT]="./xtk04 -f tk04"
    COUNT=$(( $COUNT + 1 ))
fi

# execute all interactive commands set up by previous stanzas.
for index in $INDEX_LIST ; do
    pushd ${DIRECTORY[$index]}
    echo "${COMMAND[$index]}"
    ${COMMAND[$index]} 2> test.error
    # Look for any status codes (segfaults, plexit) from the examples themselves
    status_code=$?
    if [ "$status_code" -ne 0 ] ; then
        echo "ERROR indicated by status code = $status_code for ${COMMAND[$index]}"
	OVERALL_STATUS_CODE=$status_code
    fi
    cat test.error
    # Look for any PLPLOT ERROR messages from plwarn that do not result in an exit code.
    is_error=`grep -l 'PLPLOT ERROR' test.error`
    if [ -n "$is_error" ] ; then
        echo "ERROR indicated by 'PLPLOT ERROR' in stderr for ${COMMAND[$index]}"
	OVERALL_STATUS_CODE=1
    fi
    popd
done

if [ "OFF" = "ON" ] ; then

    cd "${SRC_EXAMPLES_DIR}"/tcl
    plserver <<EOF
plstdwin .
plxframe .plw
pack append . .plw {left expand fill}
source plgrid.tcl
proc 1 {} "plgrid .plw.plwin"
1
exit
EOF

    # Look for any status codes (segfaults or whatever)
    status_code=$?
    if [ "$status_code" -ne 0 ] ; then
        echo "ERROR indicated by status code = $status_code for plgrid example"
	OVERALL_STATUS_CODE=$status_code
    fi

    if [ "$TEST_PLTCL_STANDARD_EXAMPLES" -ne 0 ] ; then
	# Copied almost entirely from pltcl_standard_examples script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tcl
        # Drop 14 because multiple devices do not seem to work in this context.
        # Drop 31 since it produces empty plot (by design).
        # Drop 32 since it has not been propagated (for the present by
        # design) from C.
	pltcl <<EOF
source tcldemos.tcl
plsdev "xwin"
plinit
# Disable pausing.
plspause 0
0
1
2
3
4
5
6
7
8
9
10
11
12
13
#14
15
16
17
18
19
20
plspause 0
21
22
23
24
25
26
27
28
29
30
#31
#32
33
exit
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for pltcl_standard_examples."
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi

    if [ "$TEST_TCLSH_STANDARD_EXAMPLES" -ne 0 ] ; then
	# Copied almost entirely from TCLSH_standard_examples script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tcl
        # Drop 14 because multiple devices do not seem to work in this context.
        # Drop 31 since it produces empty plot (by design).
        # Drop 32 since it has not been propagated (for the present by
        # design) from C.
	tclsh <<EOF
lappend auto_path "/home/pi/plplot/install_directory/share/plplot5.13.0"
package require Pltcl
source tcldemos.tcl
plsdev "xwin"
plinit
# Disable pausing.
plspause 0
0
1
2
3
4
5
6
7
8
9
10
11
12
13
#14
15
16
17
18
19
20
plspause 0
21
22
23
24
25
26
27
28
29
30
#31
#32
33
exit
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for tclsh_standard_examples."
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi

    if [ "$TEST_PLSERVER_STANDARD_EXAMPLES" -ne 0 ] ; then
	# Copied almost entirely from plserver_standard_examples script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tk
        # Drop 14 because multiple devices do not seem to work in this context.
        # Drop 31 since it produces empty plot (by design).
        # Drop 32 since it has not been propagated (for the present by
        # design) from C.
	plserver <<EOF
source tkdemos.tcl
plw::set_pause .plw 0
0
1
2
3
4
5
6
7
8
9
10
11
12
13
#14
15
16
17
18
19
20
21
22
23
24
25
26
27
28
29
30
#31
#32
33
exit
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for plserver_standard_examples"
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi
    if [ "$TEST_WISH_STANDARD_EXAMPLES" -ne 0 ] ; then
	# Copied almost entirely from wish_standard_examples script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tk
        # Drop 14 because multiple devices do not seem to work in this context.
        # Drop 31 since it produces empty plot (by design).
        # Drop 32 since it has not been propagated (for the present by
        # design) from C.
        # Both examples 2 ("Couldn't parse color 76") and 24 ("illegal
        # number of colors in cmap0: red") error out so we comment out
        # those examples for now.  Example 20 enables plspause so
        # we comment that out for now.
	wish <<EOF
lappend auto_path "/home/pi/plplot/install_directory/share/plplot5.13.0"
package require Pltk
source tkdemos.tcl
# Note wish currently disables pauses so no special method is
# required to do that (unlike the plserver case).
0
1
#2
3
4
5
6
7
8
9
10
11
12
13
#14
15
16
17
18
19
#20
21
22
23
#24
25
26
27
28
29
30
#31
#32
33
exit
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for wish_standard_examples"
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi

    if [ "$TEST_PLSERVER_RUNALLDEMOS" -ne 0 ] ; then
	# Copied almost entirely from plserver_runAllDemos script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tk
	plserver <<EOF
source runAllDemos.tcl
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for plserver_runAllDemos."
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi

    if [ "$TEST_WISH_RUNALLDEMOS" -ne 0 ] ; then
	# Copied almost entirely from wish_runAllDemos script.
        cd ..
        # Just in case EXAMPLES_DIR is a relative path such as the default ".".
	cd "${EXAMPLES_DIR}"/tk
	wish <<EOF
lappend auto_path "/home/pi/plplot/install_directory/share/plplot5.13.0"
package require Plplotter
source runAllDemos.tcl
EOF
        # Look for any status codes (segfaults or whatever)
	status_code=$?
	if [ "$status_code" -ne 0 ] ; then
            echo "ERROR indicated by status code = $status_code for wish_runAllDemos."
	    OVERALL_STATUS_CODE=$status_code
	fi
    fi
fi

if [ "$OVERALL_STATUS_CODE" -ne 0 ] ; then
    echo "A major error occurred for one of the interactive examples"
else
    echo "All interactive tests completed without any noticeable issues"
fi
exit $OVERALL_STATUS_CODE
