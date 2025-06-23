#!/bin/bash
# set -x
#

# DISPLAY=:0 sudo -u pi /usr/bin/zenity --info --timeout=30  --no-wrap --title="HH-MAA tester" --text="To start testing Press the <span color=\"green\">GREEN</span> button" &

DISPLAY=:0 export RET=`sudo -u pi /usr/bin/zenity --list --radiolist \
            --width=900 \
            --height=300 \
            --timeout=30 \
            --title="Change to repeat test" \
            --column=" " --column="Repeat" --column="Description" \
            TRUE    1 "NORMAL operation" \
            FALSE    5 "Repeat each test 5 times" \
            FALSE 200000 "Repeat 200000; Shutdown to stop (Red button)"`

# sleep 30
echo "RET=${RET};" > /tmp/loops.r

