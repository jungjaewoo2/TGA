rfile printsys bode bode_plot damp pzmap impulse rlocus rlocus_plot
rfile ssselect esort feedback zp2ss series rlocfind


disp("This file demonstates RLaB's ability for classical control system");
disp("design by going through the design of a YAW DAMPER for a Jet Transport");
disp("aircraft.");

pause();


disp("Define Jet Transport model MACH=0.8 H=40,000ft");

A=[-.0558, -.9968, .0802, .0415;
    .598,  -.115, -.0318, 0;
   -3.05,   .388, -.4650, 0;
    0,    0.0805, 1, 0];
B=[.0729,  .0001;
   -4.75,   1.23;
    1.53,   10.63;
    0   ,   0];
C=[0, 1, 0, 0;
   0, 0, 0, 1];
D=[0, 0;
   0, 0];
states=["beta", "yaw", "roll", "phi"];
inputs=["rudder", "aileron"];
outputs=["yaw-rate", "bank-angle"];
printsys(A,B,C,D,inputs,outputs,states)

disp("These are the state space matrices for a Jet Transport during cruise flight.");
disp("The model has two inputs and two outputs.  The units are radians for beta");
disp("(sideslip angle) and phi (bank angle) and radians/sec for yaw (yaw rate) and");
disp("roll (roll rate).  The rudder and aileron deflections are in degrees.");

pause();

disp("This model has one set of eigenvalues that are lightly damped.  They");
disp("correspond to what is called the Dutch Roll Mode.  We need to design");
disp("a compensator that increases the damping of these poles.  ");
disp(" ");
disp("Open Loop Eigenvalues");
damp(A); 
pzmap(A,[],[],[]); 
pause();

disp("Our design criteria is to provide damping ratio, zeta > 0.35, with natural");
disp("frequency Wn < 1.0 rad/sec.  We want to design the compensator using");
disp("classical methods. ");
disp(" ");
disp("Let's do some open loop analysis to determine possible control strategies.");

disp("Time response (we could use STEP or IMPULSE here)");
impulse(A,B,C,D); 
pause();

disp("The time responses show that the system is indeed lightly damped.  But the");
disp("time frame is much to long.  Let's look at the response over a smaller");
disp("time frame.  Define the time vector T before invoking IMPULSE.");

T = (0:20:0.2)[:];	// Define time vector from 0 to 20 secs in steps of 0.2

pause();

% Plot responses as separate graphs
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  if(exist(_ctb2_window)) {
    _save_window = _ctb2_window;
  }
  _ctb2_window = plstart(2,2);
}
if (_rlab_config.plot_support=="gnuplot")
{
   multiplot (2,2);     
}
pltitle("Input 1 Output 1"); impulse(A,B,C[1;],D[1;],1,T); 
pltitle("Input 1 Output 2"); impulse(A,B,C[2;],D[2;],1,T); 
pltitle("Input 2 Output 1"); impulse(A,B,C[1;],D[1;],2,T); 
pltitle("Input 2 Output 2"); impulse(A,B,C[2;],D[2;],2,T); 
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  _four_window = _ctb2_window;
  _ctb2_window = _save_window;
}
pause();

disp("Look at the plot from aileron (input 2) to bank_angle (output 2).  The");
disp("aircraft is oscillating around a non-zero bank angle.  Thus the aircraft");
disp("turns in response to an aileron impulse.  This behavior will be important");
disp("later.");

pause();

disp("Typically yaw dampers are designed using yaw-rate as the sensed output");
disp("and rudder as the input.  Let's look at that frequency response.");
R = bode(A,B,C[1;],D[1;],1);
bode_plot(R);
pause();

disp("From this frequency responses we see that the rudder has much effect around");
disp("the lightly damped Dutch roll mode (at 1 rad/sec). ");

disp("To make the design easier, extract of the subsystem from rudder to yaw_rate.");
</a;b;c;d/> = ssselect(A,B,C,D,1,1); // Extract system with input 1 and output 1

disp("Let's do some designs.  The simplest compensator is a gain.  We can determine");
disp("values for this gain using the root locus");
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

disp("Oops, looks like we need positive feedback.");
R = rlocus(a,b,-c,-d);
rlocus_plot(R);
pause(); 

disp("That looks better.  So using just simple feedback we can achieve a ");
disp("damping ratio of zeta=0.45.");

disp("Now its your turn, Use RLOCFIND to select the point on the root locus");
disp("with maximum damping.");

pause();
</k;poles/> = rlocfind(a,b,-c,-d);

disp("You chose gain: "+num2str(k));
damp(esort(poles).s);

pause();

disp("Let's form the closed loop system so that we can analyze the design.");
</ac;bc;cc;dc/> = feedback(a,b,c,d,[],[],[],-k);

disp("These eigenvalues should match the ones you chose.");
disp("Closed loop eigenvalues");
damp(ac);

disp("Time response using our time vector T");
impulse(ac,bc,cc,dc,1,T);  
pause();

disp("So the response looks pretty good.  Let's close the loop on the original");
disp("model and see how the response from the aileron looks.  Feedback using");
disp("input 1 and output 1 of plant.");
disp(" ");
disp("Feedback with selection vectors assumes positive feedback");
</Ac;Bc;Cc;Dc/> = feedback(A,B,C,D,[],[],[],k,[1],[1]);
disp("Closed loop eigenvalues");
damp(Ac);

disp("Time response");
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  _save_window = _ctb2_window;
  _ctb2_window = _four_window;
}
if (_rlab_config.plot_support=="gnuplot")
{
   multiplot (2,2);     
}
pltitle("Input 1 Output 1"); impulse(Ac,Bc,Cc[1;],Dc[1;],1,T); 
pltitle("Input 1 Output 2"); impulse(Ac,Bc,Cc[2;],Dc[2;],1,T); 
pltitle("Input 2 Output 1"); impulse(Ac,Bc,Cc[1;],Dc[1;],2,T); 
pltitle("Input 2 Output 2"); impulse(Ac,Bc,Cc[2;],Dc[2;],2,T); 
if (_rlab_config.plot_support=="gnuplot")
{
   nomultiplot ();     
}

if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  _ctb2_window = _save_window;
}
pause();


disp("Look at the plot from aileron (input 2) to bank_angle (output 2).  When");
disp("we move the aileron the system no longer continues to bank like a normal");
disp("aircraft.  We have over-stabilized the spiral mode.  The spiral mode is");
disp("typically a very slow mode and allows the aircraft to bank and turn without");
disp("constant aileron input.  Pilots are used to this behavior and will not like");
disp("our design if it doesn't allow them to fly normally.  Our design has moved");
disp("the spiral mode so that it has a faster frequency.");

pause();

disp("What we need to do is make sure the spiral mode doesn't move farther into");
disp("the left half plane when we close the loop.  One way flight control designers");
disp("have fixed this problem is to use a washout filter, i.e.");
disp("             Ks");
disp("   H(s) = --------");
disp("           (s + a)");
disp(" ");
disp("Choosing a = 0.333 for a time constant of 3 seconds, form the washout");
</aw;bw;cw;dw/> = zp2ss([0],[-.333],1);

disp("Connect the washout in series with our design model");
</a;b;c;d/> = series(a,b,c,d,aw,bw,cw,dw);

disp("Do another root locus");
R = rlocus(a,b,-c,-d);
rlocus_plot(R);
pause(); 

disp("Now the maximum damping is zeta = 0.25  Using RLOCFIND chose the gain for");
disp("maximum damping:");

pause();

</k;poles/> = rlocfind(a,b,-c,-d);
disp("You chose gain: "+num2str(k));
damp(esort(poles).s);

disp("Look at the closed loop response");
</ac;bc;cc;dc/> = feedback(a,b,c,d,[],[],[],-k);
impulse(ac,bc,cc,dc,1,T);  
pause();

disp("Now form the controller (washout + gain) ");
</aw;bw;cw;dw/> = series(aw,bw,cw,dw,[],[],[],k);

disp("Close the loop with the original model");
</Ac;Bc;Cc;Dc/> = feedback(A,B,C,D,aw,bw,cw,dw,[1],[1]);

disp("Final closed-loop response");
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  _save_window = _ctb2_window;
  _ctb2_window = _four_window;
}
if (_rlab_config.plot_support=="gnuplot")
{
   multiplot (2,2);     
}
pltitle("Input 1 Output 1"); impulse(Ac,Bc,Cc[1;],Dc[1;],1,T); 
pltitle("Input 1 Output 2"); impulse(Ac,Bc,Cc[2;],Dc[2;],1,T); 
pltitle("Input 2 Output 1"); impulse(Ac,Bc,Cc[1;],Dc[1;],2,T); 
pltitle("Input 2 Output 2"); impulse(Ac,Bc,Cc[2;],Dc[2;],2,T); 
if (_rlab_config.plot_support=="gnuplot")
{
   nomultiplot ();     
}
if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
{
  _ctb2_window = _save_window;
  plwin(_ctb2_window);
}
disp("Although we didn't quite meet the criteria, our design increased the damping");
disp("of the system substantially and the design does allow the pilots to fly the");
disp("aircraft normally.");



