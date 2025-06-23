//
// Simple contour plotting demonstration
//

NX = 20;
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

//
// we are running in the circles and away
//
t = [0:5:1/32]';
v = [0.5,0.1,-1];
r = 1.5;
f = 2;
xyz = v .* t + r .* [cos(2*pi*f .* t), sin(2*pi*f .* t), 0];

// Now create some plots
plwins (1,,[600,480]);

plxlabel("X-Axis");
plylabel("Y-Axis");

plimits(-3,3,-3,3,-1,1);
plgrid ("", "");
_data_3d = <<x=xx;y=yy;z=zz>>;
pltitle ("Single Surface Plot Demonstration: opt='surf'");
plsplot (_data_3d, "surf");
pause();
pltitle ("Single Surface Plot Demonstration: opt='mesh'");
plsplot (_data_3d, "mesh");
pause();

pltitle ("Two Surface and Line Plots Demonstration");
_data_3d = <<>>;
_data_3d.a1 = <<x=xx;y=yy;z=zz>>;
_data_3d.a2 = <<x=x1;y=y1;z=z1>>;
_data_3d.a3 = <<x=xyz[;1];y=xyz[;2];z=xyz[;3]>>;
plsplot (_data_3d, ["mesh","mesh",""]);
pause();

plaspect(0.5);
plimits();
pltitle ("Contour Demonstration, Aspect = 0.5");
plgrid();
plcont (<< x = xx; y = yy; z = zz>>);
plaspect ();
pause();

pltitle ("Contour Demonstration");
plimits();
plgrid ("", "");
plsplot (<<x=x1;y=y1;z=z1>> );
pause();

pltitle ("Contour Demonstration, Aspect = 1.0");
plimits();
plgrid ();
plcont (<< x = x1; y = y1; z = z1>>);
pause();
plcont (_data_3d);




