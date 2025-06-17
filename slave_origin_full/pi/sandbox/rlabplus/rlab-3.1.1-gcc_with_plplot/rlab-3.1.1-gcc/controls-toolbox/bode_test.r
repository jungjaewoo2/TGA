//
// test bode plots
//

rfile bode bode_plot

// --------------------
num = [.01, .01^2, .01];
den = [.25, .01, 1, 0, 0];
R = bode(num,den);
bode_plot(R);
pause();

// --------------------
num = [10];
den = [1,.4,4,0];
R = bode(num,den);
bode_plot(R);
pause();

// --------------------
num = [1];
den = [1,1,0];
R = bode(num,den);
bode_plot(R);
pause();

</a;b;c;d/> = tf2ss(num,den,2);
R = bode(a,b,c,d);
bode_plot(R);
pause();
 
// --------------------
num = [200,100];
den = [1,60,500,0];
R = bode(num,den);
bode_plot(R);
pause();

// --------------------
w=(.1:150:.1)';
num = [10, 10;
       10,-10]; 
den = [1,10];
R = bode(num,den,w);
bode_plot(R);
pause();


