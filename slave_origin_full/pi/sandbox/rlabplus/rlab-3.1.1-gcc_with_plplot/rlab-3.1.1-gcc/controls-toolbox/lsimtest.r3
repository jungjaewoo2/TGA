//
// This routine tests the lsim routine
//

rfile lsim

num=[2,5,1];
den=[1,2,3];
t=(0.0:10.0:0.1)';

period=4.0*ones(t.nr,t.nc);
u=(rem(t,period) >= period ./ 2);

J=lsim(num,den,u,t);

