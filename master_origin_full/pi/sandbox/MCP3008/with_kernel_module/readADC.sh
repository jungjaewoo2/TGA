#!/bin/sh
#(read channels 0 and 1 every 3 seconds)
# STATUS: UNTESTED

while true; do
        for i in 0 1; do
                echo -n "adc[${i}]: "
                cat /sys/bus/iio/devices/iio:device0/in_voltage${i}_raw
        done

        echo ""
        sleep 3
done

