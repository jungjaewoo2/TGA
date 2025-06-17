//
// test root locus plot
//

rfile printsys rlocus rlocus_plot tf2ss

num = [1];
den = [1,1,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,0,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,0];
den = [1,0,1];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,4,0,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,12,0,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,9,0,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1];
den = [1,4,9,10,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,0.2,16.01];
den = [1,0.2,25.01,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();

num = [1,0.2,25.01];
den = [1,0.2,16.01,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,9,28,40,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1];
den = [1,8,32,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,1];
den = [1,5,6,0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,4];
den = [1, 2, 0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1];
den = [1, 3, 2, 0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,2];
den = [1, 2, 2];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);
pause();

num = [1,3];
den = [2, 1, 1, 2, 0];
printsys(num,den);
R = rlocus(num,den);
rlocus_plot(R);
pause();
</a;b;c;d/> = tf2ss(num,den,2);
printsys(a,b,c,d);
R = rlocus(a,b,c,d);
rlocus_plot(R);


