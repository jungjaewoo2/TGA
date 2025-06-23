disp("This demo shows the use of some of the control system design");
disp("and analysis tools available in RLaB.");

rfile bode bode_plot conv damp lqr printmat printsys 
rfile poly2str rlocus rlocus_plot step tf2ss

disp("Suppose we start with a plant description in transfer function");
disp("form:                 2");
disp("                  .2 s  +  .3 s  +  1");
disp("       H(s)  =  ---------------------------");
disp("                  2");
disp("                (s  +  .4 s  +  1) (s + .5)");
disp(" ");
disp("We enter the numerator and denominator coefficients separately");
disp("into RLaB:");

num = [.2,  .3,  1];
den1 = [1, .4, 1];
den2 = [1, .5];
printf("num  = %s\n", poly2str(num).sout);
printf("den1 = %s\n", poly2str(den1).sout);
printf("den2 = %s\n", poly2str(den2).sout);
pause();

disp("The denominator polynomial is the product of the two terms. We");
disp("use convolution to obtain the polynomial product:");

den = conv(den1,den2);
printsys(num,den);
pause();

disp("We can look at the natural frequencies and damping factors of the");
disp("plant poles:");

damp(den);
pause();

disp("A root-locus can be obtained by using RLOCUS");


R = rlocus(num,den);
rlocus_plot(R);
pause();

disp("The plant may be converted to a state space representation");
disp("      .");
disp("      x = Ax + Bu");
disp("      y = Cx + Du");
disp(" ");
disp("using the tf2ss command:");

</a;b;c;d/> = tf2ss(num,den);
printsys(a,b,c,d);
pause();

disp("For systems described in state-space or by transfer functions,");
disp("the step response is found by using the STEP command:");

pltitle("Step response");
step(a,b,c,d,1); 
pause();

disp("The frequency response is found by using the BODE command:");

B = bode(a,b,c,d,1); 
bode_plot(B);
pause();

disp("A linear quadratic regulator could be designed for this plant.");
disp("For control and state penalties:");

r = 1
q = eye(size(a))

disp("the quadratic optimal gains, the associated Riccati solution,");
disp("and the closed-loop eigenvalues are:");


LQ = lqr(a,b,q,r);
printmat(LQ.e,"gain");
printmat(LQ.k,"solution")
printmat(LQ.s,"eigenvalues");
