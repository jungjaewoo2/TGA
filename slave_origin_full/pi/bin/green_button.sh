#!/bin/bash
# set -x
#

sudo echo "18" > /sys/class/gpio/export  # pin 12
# sudo echo "in" > /sys/class/gpio/gpio18/direction

sudo echo "17" > /sys/class/gpio/export  # pin 11
# sudo echo "out" > /sys/class/gpio/gpio17/direction
# active LOW
# sudo echo "1" > /sys/class/gpio/gpio17/value

while true
do
VAL=`cat /sys/class/gpio/gpio18/value`
if [ $VAL -eq "1" ] 
then
sudo echo "1" > /sys/class/gpio/gpio17/value
sudo /home/pi/sandbox/MCP3202_with_pigpio/prodMCP3202 > /tmp/prod.dat
sleep 0.5
# sudo -u pi lxterminal --command=/home/pi/sandbox/tf_fitting_rfiles/complete/extract_Rdc_v4.r
# sudo -u pi lxterminal  --geometry=80x24 --title=MAA --command=/home/pi/sandbox/tf_fitting_rfiles/complete/extract_Rdc_v4.r
sudo -u pi lxterminal  --geometry=80x24 --title=MAA --command=/home/pi/sandbox/tf_fitting_rfiles/complete/BEMF_plus_laser.r

for j in {1..5}
do
sudo echo "1" > /sys/class/gpio/gpio17/value
sleep 1
sudo echo "0" > /sys/class/gpio/gpio17/value
done

# play "10Hz-3kHz_chirp_Ampl=0.5.wav"
# sleep 4
fi
sleep .1
done

# never executed
sudo echo "18" > /sys/class/gpio/unexport
sudo echo "17" > /sys/class/gpio/unexport

echo "Done"

