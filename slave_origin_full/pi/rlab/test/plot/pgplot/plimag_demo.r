#
# Simple contour plotting demonstration
#

NX = 40;
NY = 40;
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

plfillstyle (2);
plgrid ("", "");
plmesh (xx, yy, zz);
pause();

plaspect(0.5);
pltitle ("Contour Demonstration, Aspect = 0.5");
plgrid();
plimag (<< x = xx; y = yy; z = zz>>);
plaspect ();
pause();

pltitle ("Contour Demonstration");
plgrid ("", "");
plmesh ( x1, y1, z1 );
pause();

pltitle ("Contour Demonstration, Aspect = 1.0");
plgrid ();
plimag (<< x = x1; y = y1; z = z1>>);
