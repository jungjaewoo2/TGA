printf("******************************\n");
printf("FH-MAA software ver. 2.1\n");
DATE=time2dstr(seconds(),"%Y-%m-%d");
TIME=time2dstr(seconds(),"%H:%M:%S");
op=reads("/var/www/html/Operator.txt");
_IP=reads("|hostname -I");
if(!isempty(_IP)){myipaddr=strsplt(_IP," ");}
printf("Today is: \t%s\n",DATE);
printf("Time is: \t%s\n",TIME);
printf("Current operator name is:\t%s\n",op);
printf("IP address is: \t%s\n",myipaddr[1]);
printf("******************************\n\n");

while(1){

	#rfile stepperTest10Twang_debug.r
break;
}
