// RLaB diary file: Extract_local_IP_address.2.diary.r. Opened Sun Mar 31 06:55:36 2019

kiwiIP=reads("|ping -n -c 1 kiwisdr.local|head -1 ")
PING kiwisdr.local (169.254.7.195) 56(84) bytes of data.t=strsplt(kiwiIP,"(")
PING kiwisdr.local   169.254.7.195) 56    84) bytes of data.   
tt=strsplt(t[2],")")
169.254.7.195   56            
ip=tt[1]
169.254.7.195length(strsplt(ip,"."))
           4  
# need to add error checking


# try another way
kiwiIP=reads("|ping -n -c 1 kiwisdr.local|head -1 ")
PING kiwisdr.local (169.254.7.195) 56(84) bytes of data.t=strsplt(kiwiIP," ")
PING             kiwisdr.local    (169.254.7.195)  56(84)           bytes            of               

data.            
tt=t[3]
(169.254.7.195)tt=strsplt(t[3])
(  1  6  9  .  2  5  4  .  7  .  1  9  5  )  
ttt=sum(tt[2:length(tt)-1])
169.254.7.195ip=sum(tt[2:length(tt)-1])
169.254.7.195length(strsplt(ip,"."))
           4  
class(ip)
stringdiary();
