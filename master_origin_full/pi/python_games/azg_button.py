import time
import RPi.GPIO as GPIO
from subprocess import call

switch1 = 18

GPIO.setmode(GPIO.BCM) # Use BCM GPIO numbers

GPIO.setup(switch1, GPIO.IN)

print "start"

while True:
	    if not GPIO.input(switch1):
		        print "Button 1 pressed"
			call(["play", "10Hz-3kHz_chirp_Ampl=0.5.wav"])
        		time.sleep(0.5)

    else:
	            pass
