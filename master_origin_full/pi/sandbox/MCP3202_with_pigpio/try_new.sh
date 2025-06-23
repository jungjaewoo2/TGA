#!/bin/sh
USAGE="Usage: $0 adc.dat"
#

if [ $# != 1 ]
then
	echo "${USAGE}" 
	exit 4
fi

FN=$1
(awk '$2 != 0 {print $1*4e-5, $2/4096*3.3 }' < $FN  |tail  -n +2 |head -n -2 ;echo;\
 awk '$5 != 0 {print $1*4e-5, $5/4096*3.3 }' < $FN  |tail  -n +2 |head -n -2 ;echo;\
 awk '$4 != 0 {print $1*4e-5, ($2 - $4)/4096*3.3 }' < $FN  |tail  -n +2 |head -n -2 ;echo;\
 awk '$4 != 0 {print $1*4e-5, $4/4096*3.3 }' < $FN  |tail  -n +2 |head -n -2) |ap
