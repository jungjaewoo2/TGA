//----------------------------------------------------------------------
//
// rlocfind
//
// syntax: </k;poles/> = rlocfind(A,B,C,D)
//         </k;poles/> = rlocfind(A,B,C,D,P)
//         </k;poles/> = rlocfind(NUM,DEN)
//         </k;poles/> = rlocfind(NUM,DEN,P)
//
// Find the root locus gains for a given set of roots.
// </K;POLES/> = rlocfind(A,B,C,D) puts up a crosshair cursor in the 
// graphics window which is used to select a pole location on an 
// existing root locus.  The root locus gain associated with this 
// point is returned in K and all the system poles for this gain are
// returned in POLES.  To use this command, the root locus for the 
// SISO state-space system (A,B,C,D) must be present in the graphics
// window.  If the system is MIMO, an error message is produced.  
// rlocfind works with both continuous and discrete linear systems.
//
// </K;POLES/> = rlocfind(NUM,DEN) is used to select a point on the 
// root locus of the polynomial system G = NUM/DEN whre NUM and DEN 
// are polynomials in descending powers of s or z.
//
// When invoked with an additional right hand argument,
// 	</K;POLES/> = rlocfind(A,B,C,D,P)
//	</K,POLES/> = rlocfind(NUM,DEN,P)
// returns a vector K of gains and the matrix of associated poles, 
// POLES. The vector K contains one element for each desired root 
// location in P.  The matrix POLES contains one column for each 
// root location in P and (length(DEN)-1) or length(A) rows.
//
// See also: rlocus, rlocus_plot
//
//----------------------------------------------------------------------
require abcdchk ss2tf tfchk roots
static (ginputx)

rlocfind = function (a,b,c,d,p)
{
  global (eps,ginput,_rlab_config,_rloc_window);

  if (nargs < 2 && nargs > 5) { error("Wrong number of arguments"); }

  // Determine which syntax is being used 
  if (nargs==2) {
    // Transfer function without poles
    </den; num/> = tfchk(a,b);
    ny = num.nr;
    nn = num.nc;
    if (ny!=1) { error("RLOCFIND must be used with SISO systems."); }
    // Get one point
    if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
    {
      if(!exist(_rloc_window)) {
        error("No root-locus plot available.");
      }
      printf("Select a point in the root locus window %i\n",_rloc_window);
      p  = ginputx(_rloc_window);
    }
    if (_rlab_config.plot_support=="gnuplot")
    {
      printf("Select a point in the root locus window %i\n",0);
      p  = ginputx(0);    
    }
    re = p.x;
    im = p.y;
    p  = re + sqrt(-1)*im;

  } else { if (nargs==3) {
    // Transfer function with poles
    </den;num/> = tfchk(a,b);
    ny = num.nr;
    nn = num.nc;
    if (ny!=1) { error("RLOCFIND must be used with SISO systems."); }
    p = c;

  } else { if (nargs==4) {
    // State space system without poles
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error(msg); }
    ny = d.nr;
    nu = d.nc;
    if (ny*nu!=1) { error("RLOCFIND must be used with SISO systems."); }
    </den;num/> = ss2tf(a,b,c,d);
    
    // Get one point
    if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
    {
      if(!exist(_rloc_window)) {
        error("No root-locus plot available.");
      }
      printf("Select a point in the root locus window %i\n",_rloc_window);
      p  = ginputx(_rloc_window);
    }
    if (_rlab_config.plot_support=="gnuplot")
    {
      printf("Select a point in the root locus window %i\n",0);
      p  = ginputx(0);    
    }
    
    re = p.x;
    im = p.y;
    p  = re + sqrt(-1)*im;

  } else {			
    // State space system with poles
    msg = abcdchk(a,b,c,d);
    if (msg!="") { error(msg); }
    ny = d.nr;
    nu = d.nc;
    if (ny*nu!=1) { error("rlocfind must be used with SISO systems."); }
    </den;num/> = ss2tf(a,b,c,d);

  }}}

  // Use root locus magnitude rule to determine gain value.  Assume the
  // root locus is based on negative feedback.
  z  = roots(num); 
  e  = roots(den);
  ny = num.nr;
  ns = num.nc;
  for (i in 1:ny ) {
    tfgain[i] = num[i;min(find(abs(num[i;])>eps))];
  }
  p = p[:];	// Make sure desired roots are a column

  np=length(p); ne = length(e); nz = length(z);
  k = prod(abs(ones(ne,1)*p.'-e*ones(1,np))) ./ prod(abs(ones(nz,1)*p'-z*ones(1,np))) ./ tfgain;
  k = abs(k);

  // Determine all the poles for each gain.
  for (i in 1:length(k)) {
    if (!finite(k[i])) {
      poles[;i] = roots(num);
    } else {
      poles[;i] = roots(den+k[i]*num);
    }
  }

  // If selecting points from root locus, plot all the roots
  #if (nargs==2||nargs==4) {
  #  status = inquire('hold');
  #  holdon();
  #  plot([real(poles),imag(poles)])
  #  if (!status) { holdoff(); }
  #  selected_point = p;	// Feedback to user
  #}

  return <<k=k; poles=poles>>;
};

ginputx = function(window)
{
  global(ginput)
  
  // read a point on the graphics window
  // ask user to click on the point

  if (exist(ginput)) {
     return ginput(window);
  } else {
     printf("Unfortunately, your plotting package does not support mouse-click input\n");
     printf("  Please enter the Real coordinate of the point: ");
     x = getline("stdin").[1];
     printf("  Please enter the Imag coordinate of the point: ");
     y = getline("stdin").[1];
     printf("--> Real=%.2g  Imag=%.2g\n",x,y);
     return <<x=x;y=y>>;     
  }
};

