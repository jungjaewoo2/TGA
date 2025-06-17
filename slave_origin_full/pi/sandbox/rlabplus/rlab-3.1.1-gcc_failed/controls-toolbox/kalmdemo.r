disp(" This file demonstrates RLaB's ability to do Kalman filtering.  ");
disp(" Both a steady state filter and a time varying filter are designed");
disp(" and simulated below.");

rfile cloop dcovar ddcgain destim dlsim dlqe dlsim parallel printsys ssdelete

pause();

disp(" Given the following discrete plant");
A=[ 1.1269,   -0.4940,    0.1129;
    1.0000,         0,         0;
         0,    1.0000,         0];
B = [-0.3832;
    0.5919;
    0.5191];
C=[ 1, 0, 0];
D=0;
printsys(A,B,C,D)

disp(" Design a Kalman filter to filter out gaussian noise added to the ");
disp(" input (i.e sensor noise).");

pause();

disp(" Let's do a steady state filter design first using the function ");
disp(" DLQE.  The function DLQE determines the optimal steady state ");
disp(" filter gain based on the process noise covariance Q, and the ");
disp(" sensor noise covariance, R.");
disp(" ");
disp(" For your information the dcgain of this plant is");
ddcgain(A,B,C,D)

disp(" and the covariance of the output due to process noise with ");
disp(" covariance of 1 is:");
dcovar(A,B,C,D,1).Y

disp(" Enter the process noise covariance:");
printf("Enter a number greater than zero (e.g. 1): ");
Q = getline("stdin").[1];

disp(" Enter the sensor noise covariance:");
printf("Enter a number greater than zero (e.g. 1): ");
R = getline("stdin").[1];

disp(" Now design the steady state filter gain");

Lss = dlqe(A,B,C,Q,R).l

pause();

disp(" The steady state Kalman filter we have just designed has the ");
disp(" following update equations:");
disp("                      -         *");
disp("  time update:        x[n+1] = Ax[n] + Bu[n]");
disp(" ");
disp("                      *      -                -");
disp("  observation update: x[n] = x[n] + L(y[n] - Cx[n] - Du[n])");
disp("  ");
disp(" Since this is a steady state filter we can combine these equations");
disp(" into one state model (the Kalman filter)");
disp("  -                -");
disp("  x[n+1] = (A-ALC) x[n] + AL y[n]");
disp(" ");
disp("  *                - ");          
disp("  y[n]   = (C-CLC) x[n] + CL y[n]");
disp("  ");
disp(" The function DESTIM produces this model (but also includes the ");
disp(" estimated states as outputs):");
</af;bf;cf;df/> = destim(A,B,C,D,Lss,[1],[1]);

disp(" Remove est. state outputs");
</af;bf;cf;df/> = ssdelete(af,bf,cf,df,[],[2:4]);
printsys(af,bf,cf,df,["u","y"],["y*"],["x1", "x2", "x3"]);

pause();

disp(" Let's see how it works by generating some data and comparing the ");
disp(" filtered response with the true response");
disp("     noise         noise");
disp("       |             |");
disp("       V             V");
disp(" u -+->O-->[Plant]-->O--> y --O->[filter]-->y*");
disp("    |                               |");
disp("    +-------------------------------+");
disp(" ");
disp(" ");
disp(" To simulate the system above, we could generate the response of ");
disp(" each part separately or we could generate both together.  To ");
disp(" simulate each seperately we would use DLSIM with the plant first");
disp(" and then the filter.  Below we illustrate how to simulate both");
disp(" together.");
disp(" ");
disp(" First, build one plant that contains the original plant and the ");
disp(" filter connected as shown in the diagram above.  Use PARALLEL and");
disp(" CLOOP.   Before connecting, add a process noise input to the plant");
disp(" that is the same as u and a sensor noise input.");
a = A; b=[B,B,zeros(3,1)]; c = C; d=[D,D,1];

disp(" Now connect the original plant input (1) with the known input (1)");
disp(" of the filter using PARALLEL.");
</at;bt;ct;dt/> = parallel(a,b,c,d,af,bf,cf,df,[1],[1],[],[]);

disp(" Put the plant output (1) into the filter sensor input (4)");
</at;bt;ct;dt/> = cloop(at,bt,ct,dt,[1],[4]);

pause();

disp(" The complete system model now has 4 inputs ");
disp("    (plant input,process noise,sensor noise,sensor input),");
disp(" 2 outputs");
disp("    (plant output with noise,filtered output),");
disp(" and 6 states.");
disp(" ");
disp(" Generate a sinusoidal input vector (known)");
t = [0:100]';
u = sin(t/5);

disp(" Generate process noise and sensor noise vectors");
rand("normal",0,1)
pnoise = sqrt(Q)*rand(t.nr,t.nc);
snoise = sqrt(R)*rand(t.nr,t.nc);

pause();

disp(" Now simulate using DLSIM");
</x;y/> = dlsim(at,bt,ct,dt,[u,pnoise,snoise,zeros(t.nr,t.nc)]);

disp(" Generate true response");
ytrue = y[;1]-snoise;

disp(" Plot comparison");
xlabel("No. of samples"); ylabel("Output"); plot([t,ytrue,y]);
xlabel("No. of samples"); ylabel("Error");  plot([t,ytrue-y[;2]]); 
pause();


disp(" Compute covariance of error");
err1 = ytrue-y[;1];
err1cov = sum(err1.*err1)/length(err1);
err2 = ytrue-y[;2];
err2cov = sum(err2.*err2)/length(err2);

disp(" Covariance of error before filtering");
err1cov

disp(" Covariance of error after filtering");
err2cov

pause();

disp(" Now let's form a time varying Kalman filter to perform the same");
disp(" task.  A time varying Kalman filter is able to perform well even");
disp(" when the noise covariance is not stationary.  However for this ");
disp(" demonstration, we will use stationary covariance.");
disp(" ");
disp(" The time varying kalman filter has the following update equations");
disp("                      -         *");
disp("  time update:        x[n+1] = Ax[n] + Bu[n]");
disp("                      -         * ");
disp("                      P[n+1] = AP[n]A' + B*Q*B'");
disp(" ");
disp("                              -       -       -1");
disp("  observation update: L[n] = P[n]C'(CP[n]C'+R)");
disp("                      *      -                   -");
disp("                      x[n] = x[n] + L[n](y[n] - Cx[n] - Du[n])");
disp("                      *               -");
disp("                      P[n] = (I-L[n]C)P[n]");
disp("                      *       *");
disp("                      y[n] = Cx[n] + Du[n]");
disp(" ");

pause();

disp(" Generate the noisy plant response");
y = dlsim(A,B,C,D,u+pnoise).y + snoise;

disp(" Now filter. We can't use DLSIM here since the system is nonlinear,");
disp(" so just implement the above equations in a loop.  Use the same inputs");
disp(" as before.");

P=B*Q*B';         // Quess for initial state covariance
x=zeros(3,1);     // Initial condition on the state
yest = zeros(t.nr,t.nc);
ycov = zeros(t.nr,t.nc); 
for (i in 1:length(t))
{
  yest[i] = C*x + D*u[i];
  ycov[i] = C*P*C';

  % Time update
  x = A*x + B*u[i];
  P = A*P*A' + B*Q*B';

  % Observation update
  L = P*C'/(C*P*C'+R);
  x = x + L*(y[i] - C*x - D*u[i]);
  P = (eye(3,3)-L*C)*P;
}

disp(" Compare true response with filtered response");
ylabel("Output"); plot([t,y-snoise,t,yest]);
ylabel("Error");  plot([t,yest-y+snoise]); 
pause();


disp(" The time varying filter also estimates the output covariance");
disp(" during the estimation.  Let's plot it to see if our filter reached");
disp(" steady state (as we would expect with stationary input noise).");
ylabel("Covar");
plot([t,ycov]);
pause();

disp(" From the covariance plot we see that the output covariance did ");
disp(" indeed reach a steady state in about 5 samples.  From then on,");
disp(" our time varying fiter has the same performance as the steady ");
disp(" state version.");
disp(" ");
disp(" Compute covariance of error");
err = y-snoise-yest;
errcov = sum(err.*err)/length(err)

pause();

disp("Let's compare the final Kalman gain matrices");
L
" "
Lss

disp(" So we see that they both obtain the same gain matrix in steady");
disp(" state.");

