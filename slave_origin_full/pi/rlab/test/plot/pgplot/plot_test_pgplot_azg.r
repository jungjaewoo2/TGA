#
# Test and demonstrate some of the plotting
# functionality in RLaB/PGplot interface.
#

#
# Test 2-D plots
#

plstart (1,2);

x = (0.25:4:.25)';

//-----------------------------------------------------------------//

pltitle ("Multiple Line Colors and Styles");
xlabel ("X-axis Label");
ylabel ("Y-axis Label");
plot ( [x, x, 2*x, 3*x] );

# subplot (4);
# subplot (2);
pltitle ("Multiple Point Colors and Styles");
plstyle ("point");
plot ( [x, x, 2*x, 3*x] );
pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

	plclose();
plstart (2,2);

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

plgrid ();  # Reset axis

//-----------------------------------------------------------------//

# Show off axis types
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

plgrid ();  # Reset axis


//-----------------------------------------------------------------//

# Data with different scales
y = (0.5:6:.2)';
pltitle("Plot Data with Different Scales");
plgrid ();
plot (<< [x,exp(-x*0.5)] ; [y, exp (y*0.25), exp(y*0.1)] >>);

# Histograms

rand("normal", 0, pi);
r = rand(2000,1);

plgrid ();
pltitle("Histogram, Default No. of Bins");
plhist (r);

pltitle("Histogram, 20 Bins");
plgrid ();
plhist (r, 20);

pltitle("Histogram, 30 Bins");
plgrid ();
plhist (r, 30);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

# Make 3D data

# Sombrero

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

# Now create some plots

pltitle("3-D Plot, Sombrero");
plgrid ("", "");
plfillstyle (2);
plcolor (2);
plmesh (xx, yy, zz );

# Sin - Cos surface

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
plcolor (3);
plmesh (x1, y1, z1, , , 1);

//-----------------------------------------------------------------//

# Now show 3-D viewing angles

pltitle("3-D Plot, Sombrero, Alt = 10, Az = 10");
plcolor (4);
plmesh (xx, yy, zz, 10, 10, 1);

pltitle("3-D Plot, Sombrero, Alt = 20, Az = 30");
plcolor (5);
plmesh (xx, yy, zz, 30, 20, 1);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

#
# Simple contour plotting demonstration
#

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

# Sin - Cos surface

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

# Now create some plots

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plcolor (2);
plmesh (xx, yy, zz);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plgrid();
plcont (<< x = xx; y = yy; z = zz>>);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plgrid ();
plfillstyle (2);
plcolor (3);
plmesh (x1, y1, z1);

xlabel("X-Axis");
ylabel("Y-Axis");
pltitle ("Contour Demonstration");
plcont (<< x = x1; y = y1; z = z1>>);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//
