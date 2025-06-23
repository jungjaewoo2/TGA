#!/bin/bash
# set -x
#


touch /var/www/html/Operator.txt
OP=`cat /var/www/html/Operator.txt`

/usr/bin/zenity --text "Enter operator name" --timeout 10 --entry --entry-text=$OP >/var/www/html/Operator.txt



DISPLAY=:0 export RET=`sudo -u pi /usr/bin/zenity --list --radiolist \
            --width=900 \
            --height=200 \
            --timeout=30 \
            --title="Select operation" \
            --column=" " --column=" " --column="Description" \
            TRUE    1 "NORMAL operation" \
            FALSE 7 "Special functions"`

if (( $RET == 7 )) 
then
DISPLAY=:0 export RET=`sudo -u pi /usr/bin/zenity --list --radiolist \
            --width=900 \
            --height=200 \
            --title="Select operation" \
            --column=" " --column=" " --column="Description" \
            FALSE    2 "Skip coarse test" \
            FALSE    5 "Repeat each test 5 times" \
            FALSE 200000 "Repeat 200000; Shutdown to stop (Red button)"`
fi

# sleep 30
echo "RET=${RET};" > /tmp/loops.r
# echo "/tmp/loops.r contains:"
# cat /tmp/loops.r

