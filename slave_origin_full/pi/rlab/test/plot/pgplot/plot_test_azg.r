//
// Test 2-D plots
//

# plwins ([1,2],"/XWIN",[600,480]);
plstart (1,2);

x = (0.25:4:.25)';

//-----------------------------------------------------------------//

pltitle ("Multiple Line Colors and Styles");
xlabel ("X-axis Label");
ylabel ("Y-axis Label");
plot ( [x, x, 2*x, 3*x] );

# subplot (4);
subplot (2);
pltitle ("Multiple Point Colors and Styles");
plstyle ("point");
plot ( [x, x, 2*x, 3*x] );
pause ("Type RETURN to continue");

//-----------------------------------------------------------------//


	plclose();
plstart(2,2);
# plwins ([2,2],"/XWIN",[600,480]);

pltitle("Plot 1, No Grid");
plgrid ("", "");
plot ([x,exp(-x*0.5)]);

pltitle("Plot 2, Grid Box Only");
plgrid ("bc", "bcv");
plot ([x,2*x,3*x,4*x]);

pltitle("Plot 3, Box, and Coordinates");
plgrid ("bcnt", "bcnvt");
plot ([x,-2*x,2*x]);

pltitle("Plot 4, ");
plgrid ("bcngst", "bcngvst");
plot ([x,x,cos(2*pi*x*0.25)]);
pause ("Type RETURN to continue");

plgrid ();  // Reset axis

//-----------------------------------------------------------------//

// Show off axis types
pltitle("Linear Axes");
plgrid();
plot ([x,exp(-x*0.5)]);

pltitle("Log X - Axis");
plgrid("bcngstl");
plot ([x,exp(-x*0.5)]);

pltitle("Log Y - Axis");
plgrid("bcngst", "bcngstlv");
plot ([x,exp(-x*0.5)]);

pltitle("Log  Axes");
plgrid("bcngstl", "bcngstlv");
plot ([x,exp(-x*0.5)]);
pause ("Type RETURN to continue");

plgrid ();  // Reset axis


//-----------------------------------------------------------------//

// Data with different scales
y = (0.5:6:.2)';
pltitle("Plot Data with Different Scales");
plgrid ();
plot (<< [x,exp(-x*0.5)] ; [y, exp (y*0.25), exp(y*0.1)] >>);

// Histograms

rand("normal", 0, pi);
r = rand(2000,1);

plgrid ();
pltitle("Histogram, Default No. of Bins");
plhist (r);

pltitle("Histogram, 30 Bins");
plgrid ();
plhist (r, 30);

pltitle("2 Histograms, 30 Bins");
plgrid ();
plhist ([r, 2*r], 30);
pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

// Some more histograms...

xlabel(); ylabel(); plgrid (); plstyle();
pltitle("Histogram, Default No. of Bins");
plhist (r);

xlabel(); ylabel(); plgrid (); plstyle("line-point");
pltitle("Histogram, Default No. of Bins");
plhist (r);

xlabel(); ylabel(); plgrid (); plstyle();
pltitle("Histogram, 30 Bins");
plgrid ();
plhist (r, 30);

xlabel(); ylabel(); plgrid (); plstyle("line-point");
pltitle("Histogram, 30 Bins");
plhist (r, 30);
pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

// Make 3D data

// Sombrero

NX = 20;
NY = 20;
xx = zeros (1, NX);
yy = zeros (1, NY);
for (i in 1:NX) { xx[i] = (i - NX/2)/(NX/2); }
for (i in 1:NY) { yy[i] = (i - NY/2)/(NY/2); }

zz = zeros (NX, NY);

for (i in 1:NX)
{
  for (j in 1:NY)
  {
    r = sqrt (xx[i]^2 + yy[j]^2);
    zz[i;j] = exp (-r * r) * cos (2*pi*r);
  }
}

// Now create some plots

pltitle("3-D Plot, Sombrero");
plgrid ("", "");
plmesh (<<x=xx;y=yy;z=zz>> );

// Sin - Cos surface

x1 = -3:3:.2;
y1 = -3:3:.2;
z1 = zeros (x1.n, y1.n);

for (i in 1:x1.n)
{
  for(j in 1:y1.n)
  {
    z1[i;j] = sin(y1[j]) * cos(x1[i]);
  }
}

pltitle("3-D Plot, Cos-Sine Surface");
plmesh (<<x=x1;y=y1;z=z1>>, , , 1);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

// Now show 3-D viewing angles

pltitle("3-D Plot, Sombrero, Default View, Alt = 60, Az = 45");
plmesh (<< x = xx; y = yy; z = zz>>);

pltitle("3-D Plot, Sombrero, Alt = 10, Az = 10");
plalt(10);
plaz(10);
plmesh (<< x = xx; y = yy; z = zz>>);

pltitle("3-D Plot, Sombrero, Alt = 20, Az = 30");
plalt(20);
plaz(30);
plmesh (<< x = xx; y = yy; z = zz>>);

pltitle("3-D Plot, Sombrero, Alt = 80, Az = 60");
plalt(80);
plaz(60);
plmesh (<< x = xx; y = yy; z = zz>>);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

//
// Simple contour plotting demonstration
//

NX = 20;
NY = 20;
xx = zeros (1, NX);
yy = zeros (1, NY);
for (i in 1:NX) { xx[i] = (i - NX/2)/(NX/2); }
for (i in 1:NY) { yy[i] = (i - NY/2)/(NY/2); }

zz = zeros (NX, NY);

for (i in 1:NX)
{
  for (j in 1:NY)
  {
    r = sqrt (xx[i]^2 + yy[j]^2);
    zz[i;j] = exp (-r * r) * cos (2*pi*r);
  }
}

// Sin - Cos surface

x1 = -3:3:.2;
y1 = -3:3:.2;
z1 = zeros (x1.n, y1.n);

for (i in 1:x1.n)
{
  for(j in 1:y1.n)
  {
    z1[i;j] = sin(y1[j]) * cos(x1[i]);
  }
}

// Now create some plots

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plalt();
plaz();
plmesh (<< x = xx; y = yy; z = zz>>);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plalt();
plaz();
plcont (<< x = xx; y = yy; z = zz>>);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plalt();
plaz();
plmesh (<< x = x1; y = y1; z = z1>>);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plalt();
plaz();
plcont (<< x = x1; y = y1; z = z1>>);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

// Now try some diff scale settings (log)

plgrid3();
pltitle("3-D Plot, Sombrero");
plmesh (<< x = xx; y = yy; z = zz>>);

plgrid3("bnstul", "bnstul");
pltitle("3-D Plot, Sombrero (log-scales)");
plmesh (<< x = xx+5; y = yy+5; z = zz>>);

plgrid3();
plaxis ("log", "log");
pltitle("3-D Plot, Cos-Sine Surface");
plmesh (<< x = x1; y = y1; z = z1>>);

plgrid3("bnstul", "bnstul");
plaxis ("log", "log");
pltitle("3-D Plot, Cos-Sine Surface (log-scales)");
plmesh (<< x = x1+5; y = y1+5; z = z1>>);

pause ("Type RETURN to continue");
