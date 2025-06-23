getIP = function(p){

	// get IP address of local link-local processor
	// function returns IPv4 IP string
	// p -> hostname string
	// [azg] Sun Mar 31 14:46:24 KST 2019



t=reads("|ping -n -c 1 "+p+".local|head -1 ");
t=strsplt(t," ");
tt=strsplt(t[3]);
ip=sum(tt[2:length(tt)-1]);
if( length(strsplt(ip,".")) != 4 ){stop("Function getIP(): Failed to get IP address of "+p+ " processor");}

return ip;
};
