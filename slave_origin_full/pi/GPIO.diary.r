// RLaB diary file: GPIO.diary.r. Opened Tue Aug 14 02:57:12 2018

 reads("|cat /sys/class/gpio/gpio18/value")
0system (" echo 1 > /sys/class/gpio/gpio17/value")
           0  
system (" echo 0 > /sys/class/gpio/gpio17/value")
           0  
# that works (above)
# the stuff below does not:
GREEN="/sys/class/gpio/gpio17/value";
OFF=["0\n"];
ON=["1\n"];
writem(GREEN,ON);
diary();
