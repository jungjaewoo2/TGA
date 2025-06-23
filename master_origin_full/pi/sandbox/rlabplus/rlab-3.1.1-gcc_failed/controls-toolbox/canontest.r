// test canon.r
rfile canon.r
format(8,5);

F = [  0,   2,  0,   0,  0;
     -.1,-.35, .1,  .1,.75;
       0,   0,  0,   2,  0;
      .4,  .4,-.4,-1.4,  0;
       0,-.03,  0,   0, -1];
G  = [0, 0, 0, 0, 1]';
H3 = [.5, 0, .5, 0, 0];
J  = 0;

can = canon(F,G,H3,J,"modal");
can.a
can.b
can.c
can.d
can.T

