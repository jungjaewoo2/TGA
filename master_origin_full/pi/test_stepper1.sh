# Example #1 (A move to absolute position)
echo "/1A12345R" > /dev/ttyUSB0
sleep 1
# Example #2 (Move loop with waits)
echo "/1gA10000M500A0M500G10R" > /dev/ttyUSB0
sleep 10
# End
echo "/1T" > /dev/ttyUSB0
echo ""

