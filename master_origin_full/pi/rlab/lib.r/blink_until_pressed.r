blink_until_pressed = function(){
# STATUS : SEEMS TO WORK
# [azg] Fri Aug 31 02:38:03 PDT 2018
#
# This is for GREEN button on the controller
////////////////////////////////////////////////////
 while( ! strtod(reads("|cat /sys/class/gpio/gpio18/value")) ){
 #while( strtod(reads("|cat /sys/class/gpio/gpio18/value")) ){

	 # OK# printf("\r Press GREEN button when ready to continue");
	 printf("\r Press ");
	 colors("green");
	 printf("GREEN ");
	 colors();
	 printf("button when ready to continue");
	 system (" echo 0 > /sys/class/gpio/gpio17/value");
	 sleep(0.5);
	 # This leaves LED ON after pressing button
	 system (" echo 1 > /sys/class/gpio/gpio17/value"); 
	 sleep(0.5);

 }
	 printf("\r                                                   \r");


 # JUST TURN IT OFF NOW
# system (" echo 0 > /sys/class/gpio/gpio17/value");
};
