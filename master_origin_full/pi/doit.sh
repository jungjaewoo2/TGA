#!/bin/bash
# set -x
#

echo "18" > /sys/class/gpio/export  # pin 12
echo "in" > /sys/class/gpio/gpio18/direction

echo "17" > /sys/class/gpio/export  # pin 11
echo "out" > /sys/class/gpio/gpio17/direction
# active LOW
# sudo echo "1" > /sys/class/gpio/gpio17/value

while true
do
VAL=`cat /sys/class/gpio/gpio18/value`
if [ $VAL -eq "1" ] 
then
for j in {1..5}
do
echo "1" > /sys/class/gpio/gpio17/value
sleep 1
echo "0" > /sys/class/gpio/gpio17/value
sleep 1
echo "1" > /sys/class/gpio/gpio17/value
sleep 1
echo "0" > /sys/class/gpio/gpio17/value
done

# play "10Hz-3kHz_chirp_Ampl=0.5.wav"
# sleep 4
fi
sleep .1
done

# never executed
echo "18" > /sys/class/gpio/unexport
echo "17" > /sys/class/gpio/unexport

echo "Done"

