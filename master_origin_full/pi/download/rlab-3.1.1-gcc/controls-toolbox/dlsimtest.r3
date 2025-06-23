//
// This routine tests the dlsim routine
//

rfile dlsim

num=[2,-3.4,1.5];
den=[1,-1.6,0.8];
u=rand(100,1);


J=dlsim(num,den,u);

