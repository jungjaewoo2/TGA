#!/bin/bash
# Test suite for c examples.
#
# Copyright (C) 2004  Alan W. Irwin
# Copyright (C) 2004  Andrew Ross
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

# This is called from plplot-test.sh with $cdir, $device, $dsuffix,
# $options, and possibly $verbose_test defined.

# To build the c examples before running this script do the following:
# pushd $cdir; make; popd

# Do the standard non-interactive examples.
lang="c"
export index lang
for index in 00 01 02 03 04 05 06 07 08 09 10 11 12 13 15 16 18 19 20 \
    21 22 23 24 25 26 27 28 30 31 33 ${critical_examples}; do
    if [ "$verbose_test" ] ; then
	echo "x${index}${lang}"
    fi
    if [ "$index" = "14" ] ; then
	echo "${OUTPUT_DIR}"/x${index}a${lang}%n.$dsuffix | \
	$DEBUG_CMD "$cdir"/x${index}${lang} -dev $device -o "${OUTPUT_DIR}"/x${index}${lang}%n.$dsuffix \
	    $options 2> test.error >| "${OUTPUT_DIR}"/x${index}${lang}_${dsuffix}.txt
    # Look for any status codes (segfaults, plexit) from the examples themselves.
	status_code=$?
    else
	$DEBUG_CMD "$cdir"/x${index}${lang} -dev $device -o "${OUTPUT_DIR}"/x${index}${lang}%n.$dsuffix \
	    $options 2> test.error >| "${OUTPUT_DIR}"/x${index}${lang}_${dsuffix}.txt
    # Look for any status codes (segfaults, plexit) from the examples themselves.
	status_code=$?
    fi
    cat test.error
    if [ "$status_code" -ne 0 ] ; then
	exit $status_code
    fi
  # Look for any PLPLOT ERROR messages from plwarn that do not result in an
  # exit code.
    is_error=`grep -l 'PLPLOT ERROR' test.error`
    if [ -n "$is_error" ] ; then
	exit 1
    fi
done
