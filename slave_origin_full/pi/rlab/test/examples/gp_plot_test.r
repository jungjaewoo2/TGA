#
# Test 2-D plots
#

x = (0.25:4:.25)';

//-----------------------------------------------------------------//

// [azg] // multiplot (1,2);

ptitle ("Multiple Line Colors and Styles");
xlabel ("X-axis Label");
ylabel ("Y-axis Label");
linestyle ("lines");
linenumbers ([1,3,6]);
plot ( [x, x, 2*x, 3*x] );

ptitle ("Multiple Point Colors and Styles");
linestyle (["lines", "linespoints", "points"]);
linenumbers (1:8);
plot ( [x, x, 2*x, 3*x] );

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

multiplot (1,2);

ptitle("Plot 1, Grid");
grid ();
plot ([x,exp(-x*0.5)]);

ptitle("Plot 2, No Grid");
nogrid ();
plot ([x,x,cos(2*pi*x*0.25)]);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

# Show off axis types

multiplot (2,2);

ptitle("Linear Axes");
grid();
nolog();
plot ([x,exp(-x*0.5)]);

ptitle("Log X - Axis");
semilogx ();
grid();
plot ([x,exp(-x*0.5)]);

ptitle("Log Y - Axis");
grid();
semilogy ();
plot ([x,exp(-x*0.5)]);

ptitle("Log  Axes");
loglog();
grid();
plot ([x,exp(-x*0.5)]);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

nolog ();

# Data with different scales
y = (0.5:6:.2)';
ptitle("Plot Data with Different Scales");
grid ();
plot (<< [x,exp(-x*0.5)] ; [y, exp (y*0.25), exp(y*0.1)] >>);

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

linestyle ("lines");
ptitle("3-D Plot, Sombrero");
pview ();
plmesh (<< x = xx; y = yy; z = zz>>);

pause ("Type RETURN to continue");

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

ptitle("3-D Plot, Cos-Sine Surface");
plmesh (<< x = x1; y = y1; z = z1>>);

pause ("Type RETURN to continue");

# Slanted plane

x2 = 1:10;
y2 = 1:10;
z2 = zeros (x2.n, y2.n);
for (i in 1:x2.n) { z2[i;] = i*ones(1,y2.n); }

ptitle("3-D Plot, Plane");
plmesh (<< x = x2; y = y2; z = z2>>);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

# Now show 3-D viewing angles

ptitle("3-D Plot, Sombrero, Default View");
pview ();
plmesh (<< x = xx; y = yy; z = zz>>);
pause ("Type RETURN to continue");

ptitle("3-D Plot, Sombrero, RotX = 30");
pview (30);
plmesh (<< x = xx; y = yy; z = zz>>);
pause ("Type RETURN to continue");

ptitle("3-D Plot, Sombrero, RotZ = 30");
pview (, 30);
plmesh (<< x = xx; y = yy; z = zz>>);
pause ("Type RETURN to continue");

ptitle("3-D Plot, Sombrero, Scale = 0.5");
pview (, , 0.5);
plmesh (<< x = xx; y = yy; z = zz>>);
pause ("Type RETURN to continue");

ptitle("3-D Plot, Sombrero, ScaleZ = 0.5");
pview (, , , 0.5);
plmesh (<< x = xx; y = yy; z = zz>>);
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

xlabel ("X-Axis");
ylabel ("Y-Axis");
ptitle ("Contour Demonstration");
pview ();
plmesh (<< x = xx; y = yy; z = zz>>);
pause ("Type RETURN to continue");

xlabel ("X-Axis");
ylabel ("Y-Axis");
ptitle ("Contour Demonstration");
plcont (<< x = xx; y = yy; z = zz >>);
pause ("Type RETURN to continue");

xlabel ("X-Axis");
ylabel ("Y-Axis");
ptitle ("Contour Demonstration");
plmesh (<< x = x1; y = y1; z = z1 >>);
pause ("Type RETURN to continue");

xlabel ("X-Axis");
ylabel ("Y-Axis");
ptitle ("Contour Demonstration");
plcont (<< x = x1; y = y1; z = z1>>);

pause ("Type RETURN to continue");

//-----------------------------------------------------------------//

// Now try some diff scale settings (log)

grid();
ptitle("3-D Plot, Sombrero");
plmesh (<< x = xx; y = yy; z = zz >>);
pause ("Type RETURN to continue");

grid ();
ptitle ("3-D Plot, Sombrero (log-scales)");
loglog ();
plmesh (<< x = xx+5; y = yy+5; z = zz >>);
pause ("Type RETURN to continue");

grid ();
nolog ();
ptitle ("3-D Plot, Cos-Sine Surface");
plmesh (<< x = x1; y = y1; z = z1 >>);
pause ("Type RETURN to continue");

grid ();
loglog ();
ptitle ("3-D Plot, Cos-Sine Surface (log-scales)");
plmesh (<< x = x1+5; y = y1+5; z = z1 >>);

pause ("Type RETURN to continue");

nolog ();
