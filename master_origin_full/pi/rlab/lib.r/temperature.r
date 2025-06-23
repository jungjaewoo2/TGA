temperature = function(addr){
//
// read temperature sensor with given address
// Wed Jan 23 15:01:54 KST 2019
///////////////////////////////////////

if(!exist(addr)) {addr="28-02131acde2aa";}

path="/sys/bus/w1/devices/"+addr+"/w1_slave";
temp="N/A";
if(isfile(path)){t=reads(path); 
r=t[2;1];
rsplt=(strsplt(r));
ln=length(rsplt);
temp=num2str(strtod(sum(rsplt[ln-4 : ln-1]))/100);
}
return temp;
};
