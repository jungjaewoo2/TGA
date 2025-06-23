#!/bin/bash
# set -x
#

sudo echo "7" > /sys/class/gpio/export  # pin 26
# sudo echo "in" > /sys/class/gpio/gpio7/direction


while true
# for i in {1..30}
do
	VAL=`cat /sys/class/gpio/gpio7/value`
	LONG=`echo ${VAL}${LONG} |cut -b1-5`
#	echo $LONG
	sleep 0.5
	if [ $LONG -eq "11111" ] 
	then
		sleep 1
		ssh gen 'sudo poweroff'
		ssh slave 'sudo poweroff'
		sleep 1 
		sudo poweroff 
		sleep 10
	fi
	sleep .1
done

# never executed
# echo "7" > /sys/class/gpio/unexport

echo "Done"
