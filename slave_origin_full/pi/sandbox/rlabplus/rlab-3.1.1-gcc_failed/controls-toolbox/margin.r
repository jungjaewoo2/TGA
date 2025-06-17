//----------------------------------------------------------------------
//
// margin
//
// syntax: </Gm;Pm;Wcg;Wcp/> = margin(mag,phase,w,option)
//
// Gain margin, phase margin, and crossover frequencies.
//
// </Gm;Pm;Wcg;Wcp/> = margin(MAG,PHASE,W,option) returns gain margin Gm,
// phase margin Pm, and associated frequencies Wcg and Wcp, given
// the Bode magnitude, phase, and frequency vectors MAG, PHASE,
// and W from a system.  Interpolation is performed between frequency
// points to find the correct values. 
//
// When invoked with option = "plot" (default),  
// margin(MAG,PHASE,W,option) plots the Bode plot with the gain and 
// phase margins marked with a vertical line.
//
// MARGIN works with the frequency response of both continuous and
// discrete systems.
//
// Example:
//  f = bode(a,b,c,d,1)  or  f = bode(num, den)
//  </Gm;Pm;Wcg;Wcp/> = margin(f.mag, f.phase, f.w, "plot")
//
// See also: bode, and dbode.
//
//----------------------------------------------------------------------
require diff nargchk

margin=function(mag,phase,w,option)
{
  global(eps,pi,_rlab_config,_bode_window)
  
  msg = nargchk(3,4,nargs);
  if (msg != "") { error(msg); }
  if (!exist(option)) { option = "plot"; }
  
  Gm = []; Pm = []; Wcg = []; Wcp = [];

  w = w[:]; // Make sure freq. is a column
  magdb = 20*log10(mag);
  logw  = log10(w);

  // Find the points where the phase wraps.
  // The following code is based on the algorithm in the unwrap function.

  cutoff = 200;					// Arbitrary value > 180
  m = phase.nr; n = phase.nc;			// Assume column orientation.
  p = mod(phase-360, 360*ones(size(phase)));    // Phases modulo 360.
  dp = [p[1;];diff(p)];				// Differentiate phases.
  jumps = (dp > cutoff) + (dp < -cutoff);	// Array of jumps locations
  jvec = (jumps!=0);

  // Find points where phase crosses -180 degrees

  pad = 360*ones(1,n);
  upcross = (([p;-pad] >= -180)&&([pad;p] < -180));
  downcross = (([p;pad] <= -180)&&([-pad;p] > -180));
  crossings = upcross + downcross;
  pvec = (crossings != 0);

  // Find points where magnitude crosses 0 db

  pad = ones(1,n);
  upcross   = (([magdb;-pad] >= 0) && ([ pad;magdb] < 0));
  downcross = (([magdb; pad] <= 0) && ([-pad;magdb] > 0));
  crossings = upcross + downcross;
  mvec = (crossings != 0);

  for (i in 1:n) {
    jloc = find(jvec[;i]!=0);
    nj = length(jloc);

    // Remove points where phase has wrapped from pvec
    if (!isempty(jloc)) { pvec[jloc;i] = zeros(nj,1); }

    ploc = find(pvec[;i] != 0);
    mloc = find(mvec[;i] != 0);

    // Find phase crossover points and interpolate gain in db and log(freq)
    // at each point.
    lambda = (-180-p[ploc-1;i]) ./ (p[ploc;i]-p[ploc-1;i]);
    gain = magdb[ploc-1;i] + lambda .* (magdb[ploc;i]-magdb[ploc-1;i]);
    freq = logw[ploc-1] + lambda .* (logw[ploc]-logw[ploc-1]);

    // Look for asymptotic behavior near -180 degrees.  (30 degree tolerance).
    // Linearly extrapolate gain and frequency based on first 2 or last 2 points.
    tol = 30;
    if (m>=2) {
      if (abs(p[1;i]+180)<tol) {
        // Starts near -180 degrees.
        lambda = (-180-p[1;i]) / (p[2;i]-p[1;i]);
        exgain = magdb[1;i] + lambda * (magdb[2;i]-magdb[1;i]);
        exfreq = logw[1] + lambda * (logw[2]-logw[1]);
        gain = [gain;exgain];  
        freq = [freq;exfreq];
      }
      if (abs(p[m;i]+180)<tol) {
        // Ends near -180 degrees.
        lambda = (-180-p[m-1;i]) / (p[m;i]-p[m-1;i]);
        exgain = magdb[m-1;i] + lambda * (magdb[m;i]-magdb[m-1;i]);
        exfreq = logw[m-1] + lambda * (logw[m]-logw[m-1]);
        gain = [gain;exgain];  
        freq = [freq;exfreq];
      }
    }
      
    if (isempty(gain)) {
      Gm  = [Gm,inf()]; 
      Wcg = [Wcg,nan()];
      ndx = [];
    else
      Gmargin = min(abs(gain));
      ndx = mini(abs(gain));
      Gm = [Gm,-gain[ndx]];    
      Wcg = [Wcg,freq[ndx]];
    }

    // Find gain crossover points and interpolate phase in degrees and log(freq)
    // at each point.
    lambda = -magdb[mloc-1;i] ./ (magdb[mloc;i]-magdb[mloc-1;i]);
    ph   = p[mloc-1;i] + lambda .* (p[mloc;i]-p[mloc-1;i]);
    freq = logw[mloc-1] + lambda .* (logw[mloc]-logw[mloc-1]);

    if (isempty(ph)) {
      Pm  = [Pm,inf()];   
      Wcp = [Wcp,nan()];
      ndx = [];
    else
      Pmargin = min(abs(ph+180));
      ndx = mini(abs(ph+180));
      Pm = [Pm,(ph[ndx]+180)];     
      Wcp = [Wcp,freq[ndx]];
    }

  }

  // Convert frequency back to rad/sec and gain back to magnitudes.
  nndx = find(finite(Wcg));
  Wcg[nndx] = 10 .^ Wcg[nndx];
  nndx = find(finite(Wcp)); 
  Wcp[nndx] = 10 .^ Wcp[nndx];  
  nndx = find(finite(Gm));
  Gm[nndx]  = 10 .^ (Gm[nndx]/20);

  // plot graph and show location of margins
  if (option == "plot") {
    // is the plot window available?
    if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
    {
      if (exist(_bode_window))
      {
         plwin(_bode_window);
      else 
         _bode_window = plstart(1,2);
      }
    }
    if (_rlab_config.plot_support=="gnuplot")
    {
       multiplot (2,1);     
    }

    // plot mag diagram
    leg = [];
    for (i in 1:mag.nc) { sprintf(s,"%d",i); leg=[leg;s]; }
    leg1 = [leg,"0dB","Gm"];
    plegend(leg1);
    plaxis("log");
    sprintf(s,"Bode Magnitude  (gain cross freq %.3g rad/s, Gm %.3g db)",Wcp,20*log10(Gm));
    pltitle(s);
    xlabel("Frequency (rad/sec)");
    ylabel("Gain dB");    
   
    plot(<<[w,20*log10(mag)];[w,zeros(w.n,1)];[Wcg,-20*log10(Gm);Wcg,0]>>);


    // 180 degree line
    if (min(min(phase)) < -180) { 
      phase180 = -180;
    else if (max(max(phase)) < 180) {
      phase180 = -180;
    else
      phase180 = 180;
    }}

    // plot phase diagram
    sprintf(s,"%d",phase180);
    leg1[leg1.n-1] = s;
    leg1[leg1.n]   = "Pm";
    plegend(leg1);
    plaxis("log");    
    sprintf(s,"Bode Phase  (phase cross freq %.3g rad/s, Pm %.3g deg)",Wcg,Pm);
    pltitle(s);
    xlabel("Frequency (rad/sec)");
    ylabel("Phase deg")
    plot(<<[w,phase];[w,phase180*ones(w.n,1)];[Wcp,Pm+phase180;Wcp,phase180]>>);
    plegend("default");
    plaxis();
    pltitle();
    xlabel("");
    ylabel("")
    if (_rlab_config.plot_support=="gnuplot")
    {
      nomultiplot ();     
    }
    
  }

  return <<gm=Gm;pm=Pm;wcg=Wcg;wcp=Wcp>>;
};

