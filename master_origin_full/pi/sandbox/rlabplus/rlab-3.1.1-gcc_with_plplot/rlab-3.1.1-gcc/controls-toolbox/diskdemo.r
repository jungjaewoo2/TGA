rfile checklist

disp("This file demonstates RLaB's ability for classical digital control");
disp("system design by going through the design of a computer hard disk");
disp("read/write head position controller.");
rfile printsys bode c2dm dbode dstep rlocus rlocus_plot rlocfind
rfile ddamp damp zgrid pzmap zp2tf series feedback margin

pause(); 

disp(" Using Newton's law, the simplest model for the read/write head has the ");
disp(" following differential equation:");
disp(" ");
disp("  I*theta_ddot + C*theta_dot + K*theta = Ki * i");
disp(" ");
disp("  where I is the inertia of the head assembly");
disp("        C is the viscous damping coefficient of the bearings");
disp("        K is the return spring constant");
disp("        Ki is the motor torque constant");
disp("        theta_ddot, theta_dot, and theta are the angular acceleration,");
disp("          angular velocity and position of the head");
disp("	i is the input current");
disp(" ");

pause();

disp(" Taking the laplace transform, the transfer function is:");
disp("                Ki");
disp("  H(s) = ----------------");
disp("         I s^2 + C s + K");
disp(" ");
disp(" Using the values I=.01 Kg m^2, C=.004 Nm/(rad/sed), K=10 Nm/rad,  and");
disp(" Ki=.05 Nm/rad form the transfer function description of this system");

I = .01; C = 0.004; K = 10; Ki = .05;
NUM = [Ki];
DEN = [I, C, K];
printsys(NUM,DEN,"s");

pause();

disp(" Our task is to design a digital controller that can be used to provide ");
disp(" accurate positioning of the read/write head.  We will do the design in the");
disp(" digital domain.  ");
disp(" ");
disp(" First we must discretize our plant since it is continuous.  RLaB has ");
disp(" several methods available for this discretization using the function C2DM.");
disp(" Let's compare all the methods and chose the 'best' one.  Use sample time");
disp(" Ts = 0.005  (5 ms)");

Ts = 0.005;
</mag;phase;w/> = bode(NUM,DEN);
</den;num/> = c2dm(NUM,DEN,Ts,"zoh");     </mzoh;pzoh;wzoh/> = dbode(num,den,Ts);
</den;num/> = c2dm(NUM,DEN,Ts,"foh");     </mfoh;pfoh;wfoh/> = dbode(num,den,Ts);
</den;num/> = c2dm(NUM,DEN,Ts,"tustin");  </mtus;ptus;wtus/> = dbode(num,den,Ts);
</den;num/> = c2dm(NUM,DEN,Ts,"prewarp",30); </mpre;ppre;wpre/> = dbode(num,den,Ts);
</den;num/> = c2dm(NUM,DEN,Ts,"matched"); </mmat;pmat;wmat/> = dbode(num,den,Ts);

disp(" Now plot the results as a comparison. ");
// check if a plot window is available 
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
 if (exist(_bode_window)) {
   plwin(_bode_window);
 else
   _bode_window = plstart(1,2);
 }
}
if (_rlab_config.plot_support=="gnuplot")
{
   multiplot (2,1);     
}
plaxis("log");
plegend("default");
xlabel("Frequency (rad/sec)"); 
ylabel("Gain db")
pltitle("c2d comparison plot");
Plt = <<[w,20*log10(mag)];[wzoh,20*log10(mzoh)];[wfoh,20*log10(mfoh)];...
       [wtus,20*log10(mtus)];[wpre,20*log10(mpre)];[wmat,20*log10(mmat)]>>;
Plt = checklist(Plt);
plot(Plt);
xlabel("Frequency (rad/sec)");
ylabel("Phase deg");
Plt = <<[w,phase];[wzoh,pzoh];[wfoh,pfoh];[wtus,ptus];[wpre,ppre];[wmat,pmat]>>;
Plt = checklist(Plt);
plot(Plt);
plaxis();
if (_rlab_config.plot_support=="gnuplot")
{
   nomultiplot ();     
}
pause();


disp(" Looking at these plots it seems that they are all pretty good matches to the");
disp(" continuous response.  However, the matched pole zero method gives a marginally");
disp(" better match to the continuous response.");

pause();

disp(" Now analyze the discrete system.");
disp("Discrete system")
printsys(num,den,"z")

disp(" Plot step response");
dstep(num,den); 
pause();

disp(" The system oscillates quite a bit.  This is probably due to very light");
disp(" damping.  We can check this by computing the open loop eigenvalues.");
disp(" ");
disp("Open loop discrete eigenvalues"); 
ddamp(den,Ts); 
zgrid("new");
pzmap(1,den); 
pause();

disp(" Note that the poles are very lightly damped and near the unit circle.");
disp(" We need to design a compensator that increases the damping of this system.");
disp(" ");
disp(" Let's try to design a compensator.  The simplest compensator is a simple gain.");
rlocus_plot(rlocus(num,den)); 
pause();

disp(" As shown in the root locus, the poles quickly leave the unit circle and go");
disp(" unstable.  We need to introduce some lead or a compensator with some zeros.");
disp(" Try the compenstor:        K(z + a)");
disp("                     D(z) = --------  where a < b");
disp("                             (z + b)");

pause();

disp(" Form compensator and connect in series with our plant");
disp(" Use a = -.85 and b = 0.");
j=sqrt(-1);
</denc;numc/> = zp2tf([.85 ]',[0]',1);
</den2;num2/> = series(num,den,numc,denc);

disp(" Lets see how the frequency response has changed.");
</mag;phase;w/>   = dbode(num,den,1);	// Use normalized frequency
</mag2;phase2;w/> = dbode(num2,den2,1,w);

disp(" Now plot a comparison plot.  Press any key after the plot ...");
// check if a plot window is available 
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
 if (exist(_bode_window)) {
   plwin(_bode_window);
 else
   _bode_window = plstart(1,2);
 }
}
if (_rlab_config.plot_support=="gnuplot")
{
   multiplot (2,1);     
}
plaxis("log");
plegend("default");
xlabel("Frequency (rad/sec)");
ylabel("Gain dB")
plot(<<[w,20*log10(mag)];[w,20*log10(mag2)]>>)

xlabel("Frequency (rad/sec)");
ylabel("Phase deg");
plot(<<[w,phase];[w,phase2];[w[1],-180;w[w.n],-180]>>);
plaxis();
if (_rlab_config.plot_support=="gnuplot")
{
   nomultiplot ();     
}
pause();

disp(" So our compensator has indeed added lead.");
disp(" ");
disp(" Now let's try the root locus again with our compensator");
zgrid("new");
rlocus_plot(rlocus(num2,den2)); 
pause();

disp(" This time the poles stay within the unit circle for some time.");
disp(" Now its your turn, Using RLOCFIND chose the poles with greatest damping");
disp(" ratio.  (The lines drawn by ZGRID show the damping ratios from z=0 to 1");
disp(" in steps of .1)");

pause();
</k;poles/> = rlocfind(num2,den2);

disp("You chose gain: "+num2str(k));
ddamp(poles,Ts);

disp(" Let's form the closed loop system so that we can analyze the design.");
</denc;numc/> = feedback(num2,den2,k,1);

disp(" These eigenvalues should match the ones you chose.");
disp("Closed loop eigenvalues");
ddamp(denc,Ts);

pause();


disp(" Closed loop time response");
dstep(numc,denc); 
pause();

disp(" So the response looks pretty good and settles in about 15 samples");
disp(" which is 15*Ts secs.");
disp(" ");
disp("Our disc drive will have a seek time > "+num2str(15*Ts)+" seconds.");
pause();

disp(" Let's now look at the robustness of our design.  The most common classical");
disp(" robustness criteria is the gain and phase margin.  The criteria is determined");
disp(" by forming a unity feedback system, calculating the Bode response and looking");
disp(" for the phase and gain crossovers.  MATLAB contains a function MARGIN that");
disp(" determines the phase and gain margin given the Bode response.");
disp(" ");
disp(" Form unity feedback system by connecting our design system with the gain");
disp(" we chose.  Leave the loop open so we can compute the open loop Bode response.");
</den2;num2/> = series(num2,den2,k,1);

disp(" Compute Bode response and margins");
</mag;phase;w/> = dbode(num2,den2,Ts);
disp(" Plot Bode plot with margins");
</Gm;Pm;Wcg;Wcp/> = margin(mag,phase,w,"plot"); 
pause();

disp(" Gain margin db, @ frequency, Phase margin, @ frequency");
Margins = [20*log10(Gm),Wcg,Pm,Wcp]

disp(" Our design is robust and can tolerate a 10 db gain increase and a 40 degree");
disp(" phase lag without going unstable.  By continuing this design process we may");
disp(" be able to find a compensator that will stabilize the open loop system and ");
disp(" allow us to reduce the seek time (more damping would allow us to reduce the");
disp(" settling time in the step response).");


