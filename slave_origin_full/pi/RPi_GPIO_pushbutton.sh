#!/bin/bash
# set -x
#
# test GPIO 18 as input
# 
# some errors if 
# GPIO configured previously
# (need root permission?) but works fine
# [azg] Fri Sep 29 17:34:55 -07 2017
#
# STATUS: works fine
####################################
# GPIO numbers should be from this list
 # 0, 1, 4, 7, 8, 9, 10, 11, 14, 15, 17, 18, 21, 22, 23, 24, 25
 
 # Note that the GPIO numbers that you program here refer to the pins
 # of the BCM2835 and *not* the numbers on the pin header. 
 # So, if you want to activate GPIO7 on the header you should be 
 # using GPIO4 in this script. Likewise if you want to activate GPIO0
 # on the header you should be using GPIO17 here.
####################################

echo "18" > /sys/class/gpio/export  # pin 12
echo "in" > /sys/class/gpio/gpio18/direction

echo "17" > /sys/class/gpio/export  # pin 11
echo "out" > /sys/class/gpio/gpio17/direction
# active LOW
echo "1" > /sys/class/gpio/gpio17/value

# while true
for i in {1..3}
do
	#echo -n "i=$i"
# echo 1 > /sys/class/gpio/gpio4/value
# echo 0 > /sys/class/gpio/gpio4/value
# cat /sys/class/gpio/gpio18/value
VAL=`cat /sys/class/gpio/gpio18/value`
echo -n $VAL
#if ( $( cat /sys/class/gpio/gpio18/value)  -eq 1 )  
#  then  play "10Hz-3kHz_chirp_Ampl=0.5.wav"
#fi
if [ $VAL -eq "1" ] 
then
#echo "true"
for j in {1..5}
do
echo "0" > /sys/class/gpio/gpio17/value
sleep 1
echo "1" > /sys/class/gpio/gpio17/value
sleep 1
done

# play "10Hz-3kHz_chirp_Ampl=0.5.wav"
# sleep 4
fi
sleep 1
# sleep .1
done

# never executed
echo "18" > /sys/class/gpio/unexport
echo "17" > /sys/class/gpio/unexport

echo "Done"

