#!/bin/bash
# set -x
#

touch /var/www/html/Operator.txt
OP=$(cat /var/www/html/Operator.txt)

/usr/bin/zenity --text "Enter operator name" --timeout 30 --entry --entry-text=$OP > /var/www/html/Operator.txt

DISPLAY=:0 export RET=$(
    sudo -u pi /usr/bin/zenity \
        --list \
        --radiolist \
        --width=400 \
        --height=160 \
        --timeout=30 \
        --hide-column=2 \
        --print-column=2 \
        --title="Select operation" \
        --column=" " --column=" " --column="Description" \
        TRUE 1 "NORMAL operation" \
        FALSE 7 "Special functions"
)

if ((RET == 7)); then
    DISPLAY=:0 export RET=$(
        sudo -u pi /usr/bin/zenity \
            --list \
            --radiolist \
            --width=400 \
            --height=250 \
            --hide-column=2 \
            --print-column=2 \
            --title="Select operation" \
            --column=" " --column=" " --column="Description" \
            FALSE 2 "Skip coarse test" \
            FALSE 5 "Repeat each test 5 times" \
            FALSE 200000 "Repeat 200000; Shutdown to stop (Red button)" \
            FALSE TIME "Set Time and Date" \
            FALSE DEL "Delete ALL Measured Data" \
            FALSE CALIB "Calibrate Coarse Motion"
    )
fi

# sleep 30
if [[ $RET != "TIME" && $RET != "DEL" && $RET != "CALIB" ]]; then
    echo "RET=${RET};" > /tmp/loops.r
fi

if [[ $RET == "TIME" ]]; then
    /home/pi/bin/dateandtime.sh
fi

if [[ $RET == "DEL" ]]; then

    if /usr/bin/zenity --question --text="Delete all measured data?"; then
        /bin/rm -f /var/www/html/*.csv /var/www/html/*.pdf /var/www/html/*.zip
    else
        echo "Deletion cancelled"
    fi
fi

if [[ $RET == "CALIB" ]]; then

    if /usr/bin/zenity --question --text="<span font-size=\"xx-large\" color=\"red\"><b>STOP!</b></span>\nCoarse Motion Recalibration Requires Management Approval. PROCEED?"; then
        HASH='3b4e878293a660215c22dff6fb8310369218b040  -'
        var1=$(echo $(zenity --forms --text="Authorization" \
            --add-password="Enter Management Password" 2> /dev/null) | shasum)
        if [ "$HASH" == "$var1" ]; then
            #echo PASSED
            /home/pi/bin/calib.sh
        else
            /usr/bin/zenity --error --text "WRONG PASSWORD" 2> /dev/null
        fi

    else
        echo "Calibration cancelled"
    fi
fi
