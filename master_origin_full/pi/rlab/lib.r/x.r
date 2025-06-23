#!/home/pi/bin/rlab2
# require last_revision
last_rev= "Sun Aug 11 21:24:56 KST 2024";
# last_rev=last_revision();
printf("******************************\n");
printf("TGA software last modified on "+last_rev+"\n");
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");
_IP=reads("|hostname -I");
if(!isempty(_IP)){myipaddr=strsplt(_IP," ");}
b=myipaddr[find(strindex(myipaddr, "169.254.") != 1)];
myipaddr=b[find(strindex(b, "2601:") != 1)]; 
printf("Today is: \t%s\n",DATE);
printf("Time is: \t%s\n",TIME);
printf("Current operator name is:\t%s\n",op);
printf("IP address is: \t%s\n",myipaddr[1]);
printf("******************************\n\n");

system("sudo echo 18 > /sys/class/gpio/export 2>/dev/null");  # pin 12
system("sudo echo 17 > /sys/class/gpio/export 2>/dev/null");  # pin 11
# [azg] Fri Mar 11 18:19:24 KST 2022
system("scp -q /home/pi/rlab/lib.r/specs.r  slave:/tmp "); 
# require blink_until_pressed
rfile blink_until_pressed;
while(1){
#INF LOOP comment hash to be removed Wed Jul 31 13:12:10 KST 2019
blink_until_pressed();
rfile stepper_twang_rev5.10.no_comments.r

}
