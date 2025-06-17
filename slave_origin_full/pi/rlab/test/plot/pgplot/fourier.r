//
// Demonstrate the effect of adding terms to a Fourier expansion
//

// We want to approximate a square wave...
// The Fourier Series for a square wave is a
// sum of odd harmonics - we will demonstrate.

// Start up a plot window...

//_plsori(1);
plwins ([2,2],"/XWIN",[600,480]);


// Compute the 1st term in the Fourier series and plot.

t = (0:10:.1)';
y = sin (t);

pltitle ("Fundamental Frequency");
plwid(10);
plfont(1);
plot ( [t,y] );
plptex ("Pen width = 10", 4, 0.4);
plptex ("Font = 1 (Normal)", 4, 0.1);
pause ();


// Now add the third harmonic to the fundamental, and plot.

y = sin (t) + sin (3*t)/3;

pltitle ("1st and 3rd Harmonics");
plwid(7);
plfont(4);
plot ( [t,y] );
plptex ("Pen width = 7", 4, 0.4);
plptex ("Font = 4 (Script)", 4, 0.1);
pause ();

// Now use the first, third, fifth, seventh, and ninth harmonics.

y = sin (t) + sin (3*t)/3 + sin (5*t)/5 + sin (7*t)/7 + sin (9*t)/9;
pltitle ("1st, 3rd, 5th, 7th, 9th Harmonics");
plwid(4);
plfont(3);
plot ( [t,y] );
plptex ("Pen width = 3", 4, 0.4);
plptex ("Font = 3 (Italic)", 4, 0.1);
pause ();

//
// Now create a matrix with rows that represent adding
// more and more terms to the series.
//

t = (0:3.14:.02);
y = zeros(10,length(t));
x = zeros(t);

for (k in 1:19:2)
{
  x = x + sin(k*t)/k;
  y[(k+1)/2;] = x;
}

//
// Now make a nice 3-D plot that shows the effect
// of adding more and more terms to the series.
//

plfont (2);
plwid(1);
pltitle ("Square Wave via Fourier Series");
ylabel ("No. Terms");
//plaz (140);
plmesh (<< x = t; y=1:10; z=y' >>);
plptex ("Pen width = 1", -1, 4.5);
plptex ("Font = 2 (Roman)", -1, 3.8);

