#!/bin/sh
#(read channels 0 and 1 every 1 seconds)
#
# NOTE: NOT TRULY DIFFERENTIAL! "-" In1 s limited to ~ 100 mv around VSS
#
# STATUS: UNTESTED

while true; do
        #for i in 0 1; do
        for i in 0 ; do
                #echo -n "adc[${i}]: "
                #OUT=`cat /sys/bus/iio/devices/iio:device0/in_voltage${i}_raw`
                OUT=`cat /sys/bus/iio/devices/iio:device0/in_voltage0-voltage1_raw`
                echo -n "adc[${i}-1]: $OUT \r"
        done

        #echo ""
        sleep .1
done

