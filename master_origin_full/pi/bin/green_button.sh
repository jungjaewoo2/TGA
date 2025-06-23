#!/bin/bash
# set -x
#

sudo echo "18" > /sys/class/gpio/export  # pin 12
# sudo echo "in" > /sys/class/gpio/gpio18/direction

sudo echo "17" > /sys/class/gpio/export  # pin 11
# sudo echo "out" > /sys/class/gpio/gpio17/direction
# active LOW
# sudo echo "1" > /sys/class/gpio/gpio17/value
sleep 10

# DISPLAY=:0 sudo -u pi /usr/bin/zenity --info --timeout=30  --no-wrap --title="HH-MAA tester" --text="To start testing Press the <span color=\"green\">GREEN</span> button" &

# DISPLAY=:0 sudo -u pi /home/pi/bin/sel.sh
# DISPLAY=:0 sudo -u pi /home/pi/bin/sel2.sh
# DISPLAY=:0 sudo -u pi /home/pi/bin/sel3.sh
DISPLAY=:0 sudo -u pi /home/pi/bin/sel4.sh



# DISPLAY=:0 export RET=`sudo -u pi /usr/bin/zenity --list --radiolist \
            # --width=600 \
            # --height=200 \
            # --timeout=30 \
            # --title="Change to repeat test" \
            # --column=" " --column="Repeat" --column="Description" \
            # TRUE    1 "NORMAL operation" \
            # FALSE    5 "Repeat each test 5 times" \
            # FALSE 200000 "Repeat 200000; Shutdown to stop (Red button)"` &
# 
# sleep 30
# echo "RET=$RET;" > /tmp/loops.r


while true
do
VAL=`cat /sys/class/gpio/gpio18/value`
# if [ $VAL -eq "1" ]  # continues restart  behaviour
if [ $VAL -eq "0" ]  # default behaviour
then
sudo echo "1" > /sys/class/gpio/gpio17/value
sleep 0.5
# DISPLAY=:0 sudo -u pi lxterminal --command=/home/pi/sandbox/complete/Jan_22_rev2/x.r
DISPLAY=:0 sudo -u pi lxterminal -t "TGA tester" --command=/home/pi/rlab/lib.r/x.r 2>/tmp/TGA_errors.LOG 

sudo echo "0" > /sys/class/gpio/gpio17/value

fi
sleep .1
done

# never executed
sudo echo "18" > /sys/class/gpio/unexport
sudo echo "17" > /sys/class/gpio/unexport

echo "Done"

