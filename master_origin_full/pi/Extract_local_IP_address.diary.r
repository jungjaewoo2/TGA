// RLaB diary file: Extract_local_IP_address.diary.r. Opened Sun Mar 31 06:00:23 2019

kiwiIP=reads("|ping -n -c 1 kiwisdr.local|head -1 ")
PING kiwisdr.local (169.254.7.195) 56(84) bytes of data.ip=substr(kiwiIP,[21:33])
169.254.7.195# verify
strsplt(ip,".")
169  254  7    195  
length(strsplt(ip,"."))
           4  
NOTE: This is unreliable/broken since number of digits in IP is variable
