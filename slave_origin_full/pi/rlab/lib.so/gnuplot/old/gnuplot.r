// ----------------------------------------
//     Time-stamp: "990420, 08:04, karst"
// ----------------------------------------
// $Header: /home/ians/dev/CVS/rlab/rlab2/misc/gnuplot.r,v 1.9 1999/04/30 13:55:55 ians Exp $
//-------------------------------------------------------------------//
//	Synopsis: RLaB graphics interface

//	The following commands are available in the RLaB graphical
//	interface using GnuPlot as the plotting device.
//	These commands are intended for use with gnuplot-3.7.
//	For details, see separate help file for each command.

//	autoscale ( axis, I )
//	autotics ( I )
//	cont ( datax, datay, datazz, keynames, I )
//	copyplot ( I1, I2 )
//	defplotwin ( I )
//	epsplot ( file, mode, I )
//	grid ( string, I )
//	hist ( data, I, keynames, I )
//	key ( x, y, z, I )
//	label ( label, x, y, z, labelnr, I )
//	latexplot ( file, mode, I )
//	linenumbers ( nrs, I )
//	linestyle ( string, I )
//	loglog ( I )
//	loglog2 ( I )
//	multiplot ( row, col, I )
//	mx2tics ( freq, I )
//	mxtics ( freq, I )
//	my2tics ( freq, I )
//	mytics ( freq, I )
//	mztics ( freq, I )
//	noautoscale ( axis, I )
//	nogrid ( I )
//	nokey ( I )
//	nolabel ( labelnr, I )
//	nolog ( I )
//	nolog2 ( I )
//	nomultiplot ( I )
//	nomx2tics ( I )
//	nomxtics ( I )
//	nomy2tics ( I )
//	nomytics ( I )
//	nomztics ( I )
//	noptime ( I )
//	notitle ( I )
//	nox2tics ( I )
//	noxtics ( I )
//	noy2tics ( I )
//	noytics ( I )
//	nozeroaxis ( axis, I )
//	noztics ( I )
//	origin ( x_origin, y_origin, I )
//	pclose ( I )
//	pend ()
//	pevery ( nr, I )
//	pformat ( string1, string2, I )
//	phelp ( string, I )
//	plot ( data, keynames, I )
//	plot2 ( data, keynames, I )
//	pointsize ( size, I )
//	preset ( I )
//	printplot ( printername, I )
//	psave ( file, data )
//	pset ( string, I )
//	psize ( x_size, y_size, I )
//	psplot ( file, mode, I )
//	pssave ( file, datax, datay, datazz )
//	pstring ( string, I )
//	ptime ( I )
//	pview ( string, I )
//	range ( startx, endx, starty, endy, startz, endz, I )
//	range2 ( startx, endx, starty, endy, startz, endz, I )
//	replot ( data, keynames, I )
//	replot2 ( data, keynames, I )
//	semilogx ( I )
//	semilogx2 ( I )
//	semilogy ( I )
//	semilogy2 ( I )
//	setterm ( term, I )
//	showplot ( I )
//	showpwin ()
//	splot ( datax, datay, datazz, keynames, I )
//	title ( string, I )
//	x2label ( string, I )
//	x2range ( start, end, I )
//	x2tics ( start, incr, end, I )
//	xlabel ( string, I )
//	xrange ( start, end, I )
//	xtics ( start, incr, end, I )
//	y2label ( string, I )
//	y2range ( start, end, I )
//	y2tics ( start, incr, end, I )
//	ylabel ( string, I )
//	yrange ( start, end, I )
//	ytics ( start, incr, end, I )
//	zeroaxis ( axis, I )
//	zlabel ( string, I )
//	zrange ( start, end, I )
//	ztics ( start, incr, end, I )

//-------------------------------------------------------------------//

// Set this so that others can tell we are using gnuplot if necessary.
_rlab_config.plot_support = "gnuplot";

//
// List to contain plot-object lists.
//

static (p);
p = <<>>;

// Store the status of axis lables.
static (FoS);
FoS = "x1y1";

// Used to keep track of the default plot window.
static (DefPlotWin);
DefPlotWin = 0;

// To be able to distinguish different RLaB2 processes, use the date
static (randprocnr);
randprocnr = getline("|date '+%H.%M.%S'").[1];

//
// Internal plot-related functions
//

static (plotl)              // Plot a list.
static (plotm)              // Plot a matrix.
static (splotm)             // splot (3D) a matrix.
static (plots)              // Send a string (command) to gnuplot.
static (replotl)            // Replot a list.
static (replotm)            // Replot a matrix.
static (plot_cmd)           // The guts of alot of the plotting commands.
static (pobj_Create)        // Create a plot-object
static (pobj_Reset)
static (pobj_Destroy)
static (pobj_TmpFileName)   // Create a temporary file name for plot-files.
static (pobj_WriteData)     // Write the plot-data (numbers) to a file.
static (pobj_PlotCmd)
static (pobj_RePlotCmd)
static (pobj_SPlotCmd)
static (pobj_PlotKeyName)
static (pobj_Plot)
static (pobj_SetRm)
static (pobj_Rm)
static (get_linetype)       // Get the linetype of the current line

//
// User interface to plot functionality
//

rrange = function(R0,N){
	xrange(-R0,R0,N);
	yrange(-R0,R0,N);
	};

defplotwin = function ( N )
{
    if (exist(N)) {
	DefPlotWin = N;
	else
	return DefPlotWin
    }
};

pobj_list = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    return p.[N];
};

plot = function ( data, keynames, N )
{
  FoS = "x1y1";
  plot_cmd ( data, keynames, N );
};

plot2 = function ( data, keynames, N )
{
    FoS = "x1y2";
    plot_cmd ( data, keynames, N );
};

replot = function ( data, keynames, N )
{
    FoS = "x1y1";
    replot_cmd ( data, keynames, N );
};

replot2 = function ( data, keynames, N )
{
    FoS = "x1y2";
    replot_cmd ( data, keynames, N );
};

copyplot = function ( I1, I2 )
{
  if (!exist(I1)) { error("One plot-windows must be specified!") }
  if (!exist(I2) && I1 == DefPlotWin)  { error("Can't copy window to itself!") }

  if (!exist(I2)) {
    I2 = I1;
    I1 = DefPlotWin;
  }

  if (!exist(p.[I1])) { error("plot-object does not exist!"); }

  if (exist(p.[I2])) {
    pobj_Reset (I2);
    else
    pobj_Create (I2);
  }

  fprintf (p.[I1].prog, "save '/usr/tmp/rlab2-tmpf-gnuplot-save'\n");

  // Delay (RLaB is a bit to fast)
#   for(i in 1:10000) { i = i; }
  fprintf (p.[I2].prog, "load '/usr/tmp/rlab2-tmpf-gnuplot-save'\n");

  // Delay (RLaB is a bit to fast)
#   for(i in 1:5000) { i = i; }
  system ("rm -f /usr/tmp/rlab2-tmpf-gnuplot-save");
};

//--------------------------------------------------------------------//

plot_cmd = function ( data, keynames, N )
{
  // Check the existence keynames
  if (exist(keynames)) 
  {
    if (class(keynames) == "num") 
    {
      if (max(size(keynames)) == 1) 
      {
	N = keynames;
	clear(keynames);
      else
	if (data.nr < data.nc) 
	{
	  data = data';
	}
	if (keynames.nr < keynames.nc) 
	{
	  keynames = keynames';
	}
	data = [data,keynames];
	clear(keynames);
	if (exist(N)) 
	{
	  if (class(N) == "string") 
	  {
	    keynames = N;
	    clear(N);
	  }
	}
      }
    }
  }
  
  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  
  // Check the existence of p.[N]
  if (max(size(N)) != 1) { error ("plot(): N must be a 1-by-1"); }

  if (exist(p.[N])) 
  {
    if (p.[N].mpp.nr == 0)
    {
      // No multiplots, safe to reset.
      pobj_Reset (N);
    else
      if (p.[N].mpp[3] == 1)
      {
	// Go ahead and reset, but only for the first plot.
	pobj_Reset (N);
      else
        p.[N].pcmd = "";
      }
    }
  else
    pobj_Create (N);
  }
    
  // Set the side to plot on
  p.[N].axes = FoS;
  
  // First plot this time
  p.[N].lpnr = 0;
  
  // Check for multiplot
  if (p.[N].mpp.nr != 0) 
  {
    printf (" Plot number %i/%i done!\n",p.[N].mpp[3],p.[N].mpp[1]*p.[N].mpp[2]);

    // Set up stuff for first plot in multiplot
    if (p.[N].mpp[3] == 1) 
    {
      fprintf (p.[N].prog, "set multiplot\n");
      fprintf (p.[N].prog, "set size %s,%s\n",1./p.[N].mpp[2],...
      1./p.[N].mpp[1]);
    }
    
    N = p.[N].mpp[4];    // This seems somewhat useless at this point ?
    fprintf (p.[N].prog, "set origin %s,%s\n",...
    p.[N].mpo[p.[N].mpp[3];1],p.[N].mpo[p.[N].mpp[3];2]);
    p.[N].mpp[3] = p.[N].mpp[3]+1;
  }
  
  // Time to plot
  if (!exist(data)) 
  {
    error ("Data does not exist");
  }

  if (class(data) == "num") 
  {
    plotm (data, N, 1, keynames);
    pobj_Plot (N);
    pobj_SetRm (N);
  else if (class(data) == "string") {
    plots (data, N);
    pobj_Plot (N);
  else if (class(data) == "list") {
    pobj_Rm (N);
    plotl (data, N, keynames);
    pobj_Plot (N);
    pobj_SetRm (N);	
  else
    error ("Invalid data for plot()");
  } } }
  
  // Check for multiplot
  if (p.[N].mpp.nr != 0) 
  {
    if (p.[N].mpp[3] == p.[N].mpp[1]*p.[N].mpp[2]+1) 
    {
      fprintf (p.[N].prog, "set nomultiplot\n");
      fprintf (p.[N].prog, "set size %s,%s\n",1,1);
      fprintf (p.[N].prog, "set origin %s,%s\n",0,0);
      p.[N].mpp = [];
      p.[N].mpo = [];
    }
  }
};
   
//--------------------------------------------------------------------//
//
// User interface to replot functionality
//

replot_cmd = function ( data, keynames, N )
{    
    // Check the existence keynames
    if (exist(keynames)) {
	if (class(keynames) == "num") {
	    if (max(size(keynames)) == 1) {
		N = keynames;
		clear(keynames);
		else
		if (data.nr < data.nc) {
		    data = data';
		}
		if (keynames.nr < keynames.nc) {
		    keynames = keynames';
		}
		data = [data,keynames];
		clear(keynames);
		if (exist(N)) {
		    if (class(N) == "string") {
			keynames = N;
			clear(N);
		    }
		}
	    }
	}
    }
    
    if (exist(data) && !exist(keynames)) {
	if (max(size(data)) == 1 && class(data) == "num") {
	    N = data;
	    clear(data);
	}
    }
    
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    
    // Check the existence of p.[N]
    if (max(size(N)) != 1) { error ("plot(): N must be a 1-by-1"); }
    if (!exist(p.[N])) {
	error ("plot-object does not exist!");
    }
    
    // Set the side to plot on
    p.[N].axes = FoS;

    // First plot this time
    p.[N].lpnr = 0;
    
    // Time to plot
    if (!exist(data)) {
	fprintf (p.[N].prog, "replot\n");
	return 0;
    }
    if (class(data) == "num") {
	p.[N].repnr = p.[N].repnr+1;
	replotm (data, N, 1, keynames);
	pobj_Plot (N);
	pobj_SetRm (N);
	else if (class(data) == "string") {
	    plots (data, N);
	    pobj_Plot (N);
	    fprintf (p.[N].prog, "replot\n");
	    else if (class(data) == "list") {
		p.[N].repnr = p.[N].repnr+1;
		replotl (data, N, keynames);
		pobj_Plot (N);
		pobj_SetRm (N);	
		else
		error ("Invalid data for replot()");
    }   }   }
};


//
// User interface to splot functionality
//

splot_cyl = function ( datara, datath, datazz, N, keynames )
{
	// check the presence of the arrays
	if (!exist(datara) || !exist(datath) || !exist(datazz)) {
		printf("\n ERROR: need three data matrices in order to do succesful plot!\n");
		return 0;
		}
    
	// check their dimensionality
	if ((datara.nc != 1 )&&(datara.nr<=1)) {
		printf("\n ERROR: need a column vector for radial coordinate!\n");
		return 0;
		}
	if ((datath.nc != 1 )&&(datath.nr<=1)) {
		printf("\n ERROR: need a column vector for angular coordinate!\n");
		return 0;
		}
	if (datath.nr != datazz.nr && datara.nr == datazz.nc){
		printf("\n ERROR: matrices do not mathc. Rows(Radial) == Cols(ZZ) && Rows(Theta)== Rows(ZZ)!\n");
		return 0;
	}
    
	// Set the default splot-object
	if (!exist(N)) { N = DefPlotWin; }

	// Check the existence of p.[N]
	if (max(size(N)) != 1) { error ("splot(): N must be a 1-by-1"); }
	if (exist(p.[N])) {
		pobj_Reset (N);
		else
		pobj_Create (N);
		}

   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
	// Delay (RLaB is a bit to fast)
# 	for(i in 1:1000) { i = i; }
	printf (" Plot number %i/%i done!\n",p.[N].mpp[3],p.[N].mpp[1]*p.[N].mpp[2]);			
	if (p.[N].mpp[3] == 1) {
	    fprintf (p.[N].prog, "set multiplot\n");
	    fprintf (p.[N].prog, "set size %s,%s\n",1./p.[N].mpp[2],1./p.[N].mpp[1]);
	}
	N = p.[N].mpp[4];
	fprintf (p.[N].prog, "set origin %s,%s\n",p.[N].mpo[p.[N].mpp[3];1],p.[N].mpo[p.[N].mpp[3];2]);
	p.[N].mpp[3] = p.[N].mpp[3]+1;
   }

   // Check how many splots to do
   nsplot = (datazz.nr/datath.nr)*(datazz.nc/datara.nr);
   if (nsplot == int(nsplot)) {

       // Time to plot
       if (class(datara) != "num") {
	   error ("Invalid data class for splot() in radial coordinate (must be num)");
	   else if (class(datath) != "num") {
	       error ("Invalid data class for splot() in angular coordinate (must be num)");
	       else if (class(datazz) != "num") {
		   error ("Invalid data class for splot() in z(Ra,Th)-matrix (must be num)");
		   else
		    fprintf (p.[N].prog, "set mapping cylindrical\n");
 		    fprintf (p.[N].prog, "set view 60, 30, 1, 1\n");
 		    fprintf (p.[N].prog, "set surface\n");
# 		    fprintf (p.[N].prog, "set nocontour\n");
		   pobj_Rm (N);
		   splotm_cyl (datara, datath, datazz, nsplot, N, 1, keynames);
		   pobj_Plot (N);
		   pobj_SetRm (N);
       }   }   }
       else
       error ("Invalid data for splot(), bad radial/angular data in z-matrix");
   }
   
   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
       if (p.[N].mpp[3] == p.[N].mpp[1]*p.[N].mpp[2]+1) {
	   fprintf (p.[N].prog, "set nomultiplot\n");
	   fprintf (p.[N].prog, "set size %s,%s\n",1,1);
	   fprintf (p.[N].prog, "set origin %s,%s\n",0,0);
	   p.[N].mpp = [];
	   p.[N].mpo = [];
       }
   }
};

splot = function ( datax, datay, datazz, keynames, N )
{
    // Check what has been sent here
    // One input!
    if (!exist(datay)) {
	datazz = datax;
	datax = 1:datazz.nc;
	datax = datax';
	datay = 1:datazz.nr;
	datay = datay';
	// Two inputs!
	else 
	if (!exist(datazz)) {
	    if (exist(datay)) {
		if (class(datay) == "num") {
		    N = datay;
		    else
		    keynames = datay;
		}
	    }
	    datazz = datax;
	    datax = 1:datazz.nc;
	    datax = datax';
	    datay = 1:datazz.nr;
	    datay = datay';
	    // Three inputs!
	    else 
	    if (!exist(keynames)) {
		if (class(datay) != "num") {
		    keynames = datay;
		    N = datazz;
		    datazz = datax;
		    datax = 1:datazz.nc;
		    datax = datax';
		    datay = 1:datazz.nr;
		    datay = datay';
		}
		// Four inputs!
		else 
		if (!exist(N)) {
		    if (class(keynames) == "num") {
			N = keynames;
			clear(keynames);
		    }
    }   }   }   }
    
    // Transpose x- and y-vectors and z-matrix if nessecary
    if (datax.nc != 1 ) {
	if (datax.nr != 1 ) {
	    error ("Invalid data for x-axis (not a vector)");
	    else
	    datax = datax';
	}
    }
    if (datay.nc != 1 ) {
	if (datay.nr != 1 ) {
	    error ("Invalid data for y-axis (not a vector)");
	    else
	    datay = datay';
	}
    }
    if (datax.nr != datazz.nr && datax.nr == datazz.nc){
	datazz = datazz';
    }
    
    // Set the default splot-object
    if (!exist(N)) { N = DefPlotWin; }

    // Check the existence of p.[N]
    if (max(size(N)) != 1) { error ("splot(): N must be a 1-by-1"); }
    if (exist(p.[N])) {
	pobj_Reset (N);
	else
	pobj_Create (N);
    }

   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
	// Delay (RLaB is a bit to fast)
# 	for(i in 1:1000) { i = i; }
	printf (" Plot number %i/%i done!\n",p.[N].mpp[3],p.[N].mpp[1]*p.[N].mpp[2]);			
	if (p.[N].mpp[3] == 1) {
	    fprintf (p.[N].prog, "set multiplot\n");
	    fprintf (p.[N].prog, "set size %s,%s\n",1./p.[N].mpp[2],1./p.[N].mpp[1]);
	}
	N = p.[N].mpp[4];
	fprintf (p.[N].prog, "set origin %s,%s\n",p.[N].mpo[p.[N].mpp[3];1],p.[N].mpo[p.[N].mpp[3];2]);
	p.[N].mpp[3] = p.[N].mpp[3]+1;
   }

   // Check how many splots to do
   nsplot = (datazz.nr/datax.nr)*(datazz.nc/datay.nr);
   if (nsplot == int(nsplot)) {
       // Time to plot
       if (class(datax) != "num") {
	   error ("Invalid data class for splot() in x-vector (must be num)");
	   else if (class(datay) != "num") {
	       error ("Invalid data class for splot() in y-vector (must be num)");
	       else if (class(datazz) != "num") {
		   error ("Invalid data class for splot() in x-matrix (must be num)");
		   else
# 		    fprintf (p.[N].prog, "set view 60, 30, 1, 1\n");
# 		    fprintf (p.[N].prog, "set surface\n");
# 		    fprintf (p.[N].prog, "set nocontour\n");
		   pobj_Rm (N);
		   splotm (datax, datay, datazz, nsplot, N, 1, keynames);
		   pobj_Plot (N);
		   pobj_SetRm (N);
       }   }   }
       else
       error ("Invalid data for splot(), bad x/y-range in z-matrix");
   }
   
   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
       if (p.[N].mpp[3] == p.[N].mpp[1]*p.[N].mpp[2]+1) {
	   fprintf (p.[N].prog, "set nomultiplot\n");
	   fprintf (p.[N].prog, "set size %s,%s\n",1,1);
	   fprintf (p.[N].prog, "set origin %s,%s\n",0,0);
	   p.[N].mpp = [];
	   p.[N].mpo = [];
       }
   }
};



//
// User interface to splot functionality as 2D-contour
//

cont = function ( datax, datay, datazz, keynames, N )
{
    // Check what has been sent here
    // One input!
    if (!exist(datay)) {
	datazz = datax;
	datax = 1:datazz.nc;
	datax = datax';
	datay = 1:datazz.nr;
	datay = datay';
	// Two inputs!
	else if (!exist(datazz)) {
	    if (exist(datay)) {
		if (class(datay) == "num") {
		    N = datay;
		    else
		    keynames = datay;
		}
	    }
	    datazz = datax;
	    datax = 1:datazz.nc;
	    datax = datax';
	    datay = 1:datazz.nr;
	    datay = datay';
	    // Three inputs!
	    else if (!exist(keynames)) {
		if (class(datay) != "num") {
		    keynames = datay;
		    N = datazz;
		    datazz = datax;
		    datax = 1:datazz.nc;
		    datax = datax';
		    datay = 1:datazz.nr;
		    datay = datay';
		}
		// Four inputs!
		else if (!exist(N)) {
		    if (class(keynames) == "num") {
			N = keynames;
			clear(keynames);
		    }
    }   }   }   }
    
    // Transpose x- and y-vectors and z-matrix if nessecary
    if (datax.nc != 1 ) {
	if (datax.nr != 1 ) {
	    error ("Invalid data for x-axis (not a vector)");
	    else
	    datax = datax';
	}
    }
    if (datay.nc != 1 ) {
	if (datay.nr != 1 ) {
	    error ("Invalid data for y-axis (not a vector)");
	    else
	    datay = datay';
	}
    }
    if (datax.nr != datazz.nr && datax.nr == datazz.nc){
	datazz = datazz';
    }
    
    // Set the default splot-object
    if (!exist(N)) { N = DefPlotWin; }

    // Check the existence of p.[N]
    if (max(size(N)) != 1) { error ("splot(): N must be a 1-by-1"); }
    if (exist(p.[N])) {
	pobj_Reset (N);
	else
	pobj_Create (N);
    }

   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
	// Delay (RLaB is a bit to fast)
# 	for(i in 1:1000) { i = i; }
	printf (" Plot number %i/%i done!\n",p.[N].mpp[3],p.[N].mpp[1]*p.[N].mpp[2]);
	if (p.[N].mpp[3] == 1) {
	    fprintf (p.[N].prog, "set multiplot\n");
	    fprintf (p.[N].prog, "set size %s,%s\n",1/p.[N].mpp[2],1/p.[N].mpp[1]);
	}
	N = p.[N].mpp[4];
	fprintf (p.[N].prog, "set origin %s,%s\n",p.[N].mpo[p.[N].mpp[3];1],p.[N].mpo[p.[N].mpp[3];2]);
	p.[N].mpp[3] = p.[N].mpp[3]+1;
   }

    // Check how many splots to do
    nsplot = (datazz.nr/datax.nr)*(datazz.nc/datay.nr);
    if (nsplot == int(nsplot)) {
	// Time to plot
	if (class(datax) != "num") {
	    error ("Invalid data class for splot() in x-vector (must be num)");
	    else if (class(datay) != "num") {
		error ("Invalid data class for splot() in y-vector (must be num)");
		else if (class(datazz) != "num") {
		    error ("Invalid data class for splot() in x-matrix (must be num)");
		    else
		    fprintf (p.[N].prog, "set view 0, 0, 1\n");
		    fprintf (p.[N].prog, "set nosurface\n");
		    fprintf (p.[N].prog, "set contour\n");
		    pobj_Rm (N);
		    splotm (datax, datay, datazz, nsplot, N, 1, keynames);
		    pobj_Plot (N);
		    pobj_SetRm (N);
	}   }   }
	else
	error ("Invalid data for splot(), bad x/y-range in z-matrix");
    }
    
   // Check for multiplot
   if (p.[N].mpp.nr != 0) {
	if (p.[N].mpp[3] == p.[N].mpp[1]*p.[N].mpp[2]+1) {
	    fprintf (p.[N].prog, "set nomultiplot\n");
	    fprintf (p.[N].prog, "set size %s,%s\n",1,1);
	    fprintf (p.[N].prog, "set origin %s,%s\n",0,0);
	    p.[N].mpp = [];
	    p.[N].mpo = [];
	}
   }
};



//
// User interface to plot histograms
//

hist = function ( data, I, keynames, N )
{
    // Set the default number of bins
    if (!exist(I)) { I = 10; }
    // Three inputs
    if (exist(keynames)) {
	if (class(keynames) == "num") {
	    N = keynames;
	    clear(keynames);
	}
    }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    // sort data
    data = sort(data).val;
    for (col in 1:data.nc) {
	i = 1;
	matrix.[col] = [];
	minval = min(min(data[;col]));
	maxval = max(max(data[;col]));
	stepval = (maxval-minval)/(I);
	matrix.[col][i;1] = minval + (i-1+0.5)*stepval;
	matrix.[col][i;2] = 0;
	for (j in 1:data.nr) {
	    if (data[j;col] >= matrix.[col][i;1]-0.5*stepval && data[j;col] <= matrix.[col][i;1]+0.5*stepval) {
		matrix.[col][i;2] = matrix.[col][i;2] + 1;
		else
		if (data[j;col] >= matrix.[col][i;1]+0.5*stepval) {
		    i++;
		    matrix.[col][i;1] = minval + (i-1+0.5)*stepval;
		    matrix.[col][i;2] = 0;
		}
	    }
	}
    }
    minval = min(min(data));
    maxval = max(max(data));
    stepval = (maxval-minval)/(I);
    xrange ((minval-stepval)[1],(maxval+stepval)[1],N);
    ymax = 0;
    for (i in 1:size(matrix)) {
	newymax = max(max(matrix.[i][;2]));
	if (newymax > ymax) {ymax = newymax;}
    }
    yrange (0,int(newymax+1)[1],N);
    if (exist(keynames)) {
	plot (matrix,keynames,N);
	else
	plot (matrix,N);
    }
};

//
// Send string to GNUPLOT
//

pstring = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "%s\n",string);
};



//
// User interface to save a plot ot splot as postscript hardcopy
//

psplot = function ( file, mode, N )
{
  //
  // If mode is a number, and N is specified,
  // then use mode as N.
  //
  
  if (!exist(N) && exist(mode)) {
    if (class(mode) == "num") {
      N = mode;
    }
  }

  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  if (!exist(p.[N])) { 
    error ("No existing plot to make hardcopy from"); 
  }
  
  //
  // Check the mode... use 'default enh' if none specified.
  //
  if (!exist(mode)) { mode = "color enh"; }

  // Make hardcopy file
  if (class(mode) == "num") { mode = "color enh"; }
  fprintf (p.[N].prog, "set term post %s\n", mode);
  fprintf (p.[N].prog, "set output \"%s\"\n", file);
  fprintf (p.[N].prog, "replot\n");
  
  // Reset to original term type, and replot
  fprintf (p.[N].prog, "set term %s\n", p.[N].term);
  fprintf (p.[N].prog, "replot\n");

};

postplot = psplot;

epsplot = function ( file, mode, N )
{
  //
  // If mode is a number, and N is specified,
  // then use mode as N.
  //
  if (!exist(N) && exist(mode))
  {
    if (class(mode) == "num")
    {
      N = mode;
    }
  }

  //
  // name
  //
  if (!exist(file))
  {
    file = "./tmp.eps";
  }
  if (class(file)!="string")
  {
    file = "./tmp.eps";
  }


  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  if (!exist(p.[N])) { error ("No existing plot to make hardcopy from"); }

  //
  // Check the mode... use 'default enh' if none specified.
  // set terminal postscript eps color lw 15 "Helvetica" 20
  if (!exist(mode)) { mode = "enh"; }

  // Make hardcopy file
  if (class(mode) == "num") { mode = "enh"; }

  printf ("Printing to 'eps' file ");

  //fprintf (p.[N].prog, "set term postscript eps %s  color blacktext solid\n", mode);
  //pset ("term postscript eps enh color blacktext solid", N);
  pset ("term postscript eps enh monochrome blacktext solid lw 2 \"Times-Roman\" 21", N);
  pset ("output \""+ file + "\" \n", N);

  fprintf (p.[N].prog, "replot \n");

  printf(". Done!\n");

  // Reset to original term type, and replot
  //fprintf (p.[N].prog, "set term %s\n", p.[N].term);
  //pset ("term " + p.[N].term, N);
  //fprintf (p.[N].prog, "replot\n");

};


//
// User interface to save a plot ot splot as LaTeX hardcopy
//

latexplot = function ( file, mode, N )
{
  if (!exist(N) && exist(mode)) {
    if (class(mode) == "num") {
      N = mode;
    }
  }

  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  if (!exist(p.[N])) { error ("No existing plot to make hardcopy from"); }
  
  // Make hardcopy file
  if (!exist(mode)) { mode = "default"; }
  if (class(mode) == "num") { mode = "default"; }
  fprintf (p.[N].prog, "set term latex %s\n", mode);
  fprintf (p.[N].prog, "set output \"%s\"\n", file);
  fprintf (p.[N].prog, "replot\n");

  // Reset to original term type, and replot
  fprintf (p.[N].prog, "set term %s\n", p.[N].term);
  fprintf (p.[N].prog, "replot\n");
};

//
// User interface to print a plot ot splot
//

printplot = function ( printername, N )
{
  if (exist(printername) && !exist(N)) {
    if (class(printername) == "num") {
      N = printername;
      clear(printername);
    }
  }

  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  // Check the existence of p.[N]
  if (!exist(p.[N])) { pobj_Create (N); }

  // Create tmp-file-name
  fn = tmpnam ();

  psplot (fn, "", N);
  wait(2);
  if (exist(printername)) 
  {
    system ("cat " + fn + "| lpr -P" + printername);
  else
    system ("cat " + fn + "| lpr");
  }

  //
  // Remove the temporary file
  //
  system ("rm -f " + fn);
};


//
// User interface to save plots as GNUPLOT-datafiles
//

psave = function ( file, data )
{
    // Check file-name
    if (class(file) != "string") {
	error ("Invalid string for file_name");
    }

    // Determine how many lines to save
    nplot = max([1, data.nc - 1]);

    if (nplot > data.nr) {
	printf (" Save %i columns, are you sure [y(es)/n(o)/t(ranspose) (def. no)] ? ", data.nc);
	ans = getline ("stdin");
	if (ans.[1] != "y" && ans.[1] != "t") {
	    return 0;
	}
	if (ans.[1] == "t") {
	    data = data';
	    nplot = max([1, data.nc - 1]);
	}
    }

    // Generate two column matrices for saving
    for (i in 1:nplot) {

	// Write data to file
	if (nplot == 1) {
	    pobj_WriteData (real(data), file);
	    else
	    pobj_WriteData (real(data[;1,i+1]), file);
	}
    }
};



//
// User interface to save splots as GNUPLOT-datafiles
//

pssave = function ( file, datax, datay, datazz )
{
    // Determine how many plots to draw
    if (!exist(datazz)) { 
	datazz = datax;
	datax = 1:datazz.nr;
	datax = datax';
	datay = 1:datazz.nc;
	datay = datay';
    }

    if (datax.nc != 1 ) {
	if (datax.nr != 1 ) {
	    error ("Invalid data for x-axis (not a vector)");
	    else
	    datax = datax';
	}
    }
    if (datay.nc != 1 ) {
	if (datay.nr != 1 ) {
	    error ("Invalid data for y-axis (not a vector)");
	    else
	    datay = datay';
	}
    }

    // Create file-name
    if (!exist(file)) { error ("Invalid string for file_name"); }
    if (class(file) != "string") { error ("Invalid string for file_name"); }

    nsplot = (datazz.nr/datazz.nc)/(datax.nr/datay.nr);

    // Generate three column matrices for splot program
    // Write data to new data-file
    for (i in 1:nsplot) {
	newdata = [];
	for (m in 1:datax.nr){
	    rownumber = 0;
	    for (n in 1:datay.nr){
		rownumber = rownumber+1;
		// x-value
		newdata[rownumber;1] = datax[m];
		// y-value
		newdata[rownumber;2] = datay[n];
		// z-value
		newdata[rownumber;3] = datazz[m;n+(i-1)*datax.nr];
	    }
# 	    writem (fn,real(newdata[;1,2,3]));
# 	    writem (fn,"");
#	    // Add an empty line
# 	    close (fn);
# 	    system(" echo '' >> " + fn);
# 	    open (fn,"append");
	    a=real(newdata[;1,2,3]);
	    for(row in 1:size(a)[1]) {
		for(column in 1:size(a)[2]) {
		    fprintf(fn, a[row;column]);
		    fprintf(fn, "\t");
		}
		fprintf(fn, "\n");
	    }
	    fprintf(fn, "\n");
	}
	close (fn);
    }
};



//--------------------------------------------------------------------//
//
// Make multiplots
//

multiplot = function ( row, col, N )
{
  // Check arguments for usage.

  // Only 1 argument specified.... go into multiplot mode.
  if (exist(row) && !exist(col)) { N = row; }

  // Possibly multiple arguments, set default window.
  if (!exist(N)) { N = DefPlotWin; }

  // Check the existence of p.[N]
  if (!exist(p.[N])) { pobj_Create (N); }

  // Check if to set multiplot or to do a multiple plot.
  if (!exist(col)) 
  {
    // Go into multiplot mode, and let the user handle things.
    fprintf (p.[N].prog, "set multiplot\n");
  else 
    // Manange the plot window automatically.
    p.[N].mpp = [row,col,1,N];
    p.[N].mpo = [];

    // printf (" Make %i plots!\n",p.[N].mpp[1]*p.[N].mpp[2]);
    // Build a matrix containing the size and origin information
    // for doing the multiplot.
    k = 1;
    for (i in 1:row) 
    {
      for (j in 1:col) 
      {
	p.[N].mpo[k;1] = (j-1)/col;
	p.[N].mpo[k;2] = (row-i)/row;
	k = k+1;
      }
    }
    p.[N].mpo;
  }
};



//--------------------------------------------------------------------//
//
// End multiplot mode
//

nomultiplot = function ( N )
{
   // Set the default plot-object
   if (!exist(N)) { N = DefPlotWin; }
   // Check the existence of p.[N]
   if (!exist(p.[N])) { pobj_Create (N); }
   fprintf (p.[N].prog, "set size %s,%s\n",1,1);
   fprintf (p.[N].prog, "set origin %s,%s\n",0,0);
   fprintf (p.[N].prog, "set nomultiplot\n");
   p.[N].mpp = [];
   p.[N].mpo = [];
};


//
// User interface to reset a GNUPLOT-window
//
pset = function ( string, N )
{
  // Set the default plot-object
  if (!exist(N))
  {
    N = DefPlotWin;
  }

  // Check the existence of p.[N]
  if (!exist(p.[N]))
  {
    pobj_Create (N);
  }
  if (!exist(string))
  {
    error("String missing!");
  }
  if (class(string) != "string")
  {
    error("String must be a string!");
  }

  fprintf (p.[N].prog, "set %s\n", string);
};

preset = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    if (N == "all") {
	plwins = members(p);
	for (i in 1:plwins.nc) {
	    N = eval(plwins[i]);
	    if (exist(p.[N])) {
		fprintf (p.[N].prog, "reset\n");
	    }
	}
	N = "all";
	else
	if (!exist(p.[N])) { return -1; }
	fprintf (p.[N].prog, "reset\n");
    }
    return N;
};


//
// User interface to close a GNUPLOT-window
//

pclose = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    if (N == "all") {
	plwins = members(p);
	for (i in 1:plwins.nc) {
	    N = eval(plwins[i]);
	    if (exist(p.[N])) {
		close (p.[N].rmf);
		close (p.[N].prog);
		pobj_Destroy (N);
	    }
	}
	N = "all";
	else
	if (!exist(p.[N])) { return -1; }
	close (p.[N].rmf);
	close (p.[N].prog);
	pobj_Destroy (N);
    }
    return N;
};

pend = function ()
{
    plwins = members(p);
    for (i in 1:plwins.nc) {
	N = eval(plwins[i]);
	if (exist(p.[N])) {
	    close (p.[N].rmf);
	    close (p.[N].prog);
	    pobj_Destroy (N);
	}
    }
    DefPlotWin = 0;
};



//
// Set the terminal type for the I-th GNUPLOT process
//

setterm = function ( term, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    p.[N].term = term;

    // Now send the "set term" command to GNUPLOT
    fprintf (p.[N].prog, "set term %s\n", term);
};

showterm = function ( N )
{
  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  
  // Check the existence of p.[N]
  if (!exist(p.[N])) { pobj_Create (N); }
  
  return p.[N].term;
};

showoutput = function ( N )
{
  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }
  
  // Check the existence of p.[N]
  if (!exist(p.[N])) { pobj_Create (N); }
  
  return p.[N].output;
};

//
// Print out an element of the plot list
//

showplot = function ( N )
{
  // Set the default plot-object
  if (!exist(N)) { N = DefPlotWin; }

  for (i in N)
  {
    if (!exist(p.[i])) 
    {
      pobj_Create (i);
    }
    
    // Handle members() output if necessary
    if (class (i) == "string")
    {
      printf ("\tPlot List %s\n", i);
    else
      printf ("\tPlot List %i\n", i);
    }
    printf ("\t\tTerm:\t\t\t%s\n", p.[i].term);
    printf ("\t\tOutput:\t\t\t%s\n", p.[i].output);
    printf ("\t\tTmp Files:\t\t%s\n\n", p.[i].files);
  }
};


//
// Return the member names of the plot-list.
//

showpwin = function ()
{
  return members (p);
};


//
// Here comes a bunch of functions that really is not nessecary.
// They could be done with string("set key") etc, but this is more convinient.
//

phelp = function ( string, N )
{
    // Set the default plot-object
    if (exist(string) && !exist(N)) {
	if (class(string) == "num") {
	    N = string;
	    string = "";
	}
    }
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(string)) { string = ""; }
    fprintf (p.[N].prog, "help %s\n",string);
};

xrange = function ( start, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(end)) {
	fprintf (p.[N].prog, "set xrange [%s:] writeback\n",start);
    }
    if (!exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set xrange [:%s] writeback\n",end);
    }
    if (exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set xrange [%s:%s]\n",start,end);
    }
};

yrange = function ( start, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(end)) {
	fprintf (p.[N].prog, "set yrange [%s:] writeback\n",start);
    }
    if (!exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set yrange [:%s] writeback\n",end);
    }
    if (exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set yrange [%s:%s]\n",start,end);
    }
};

zrange = function ( start, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(end)) {
	fprintf (p.[N].prog, "set zrange [%s:] writeback\n",start);
    }
    if (!exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set zrange [:%s] writeback\n",end);
    }
    if (exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set zrange [%s:%s]\n",start,end);
    }
};


range = function ( startx, endx, starty, endy, startz, endz, N )
{
    // Check what has been sent
    if (exist(startz) && !exist(endz)) { N = startz; }
    if (exist(starty) && !exist(endy)) { N = starty; }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set xrange [%s:%s]\n",startx,endx);
    if (exist(endy)) {
	fprintf (p.[N].prog, "set yrange [%s:%s]\n",starty,endy);
    }
    if (exist(endz)) {
	fprintf (p.[N].prog, "set zrange [%s:%s]\n",startz,endz);
    }
};

x2range = function ( start, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(end)) {
	fprintf (p.[N].prog, "set x2range [%s:] writeback\n",start);
    }
    if (!exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set x2range [:%s] writeback\n",end);
    }
    if (exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set x2range [%s:%s]\n",start,end);
    }
};

y2range = function ( start, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (!exist(end)) {
	fprintf (p.[N].prog, "set y2range [%s:] writeback\n",start);
    }
    if (!exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set y2range [:%s] writeback\n",end);
    }
    if (exist(start) && exist(end)) {
	fprintf (p.[N].prog, "set y2range [%s:%s]\n",start,end);
    }
};

range2 = function ( startx, endx, starty, endy, startz, endz, N )
{
    // Check what has been sent
    if (exist(startz) && !exist(endz)) { N = startz; }
    if (exist(starty) && !exist(endy)) { N = starty; }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set x2range [%s:%s]\n",startx,endx);
    if (exist(endy)) {
	fprintf (p.[N].prog, "set y2range [%s:%s]\n",starty,endy);
    }
    if (exist(endz)) {
	fprintf (p.[N].prog, "set zrange [%s:%s]\n",startz,endz);
    }
};

pview = function( rotx, rotz, scale, scalez, N )
{
  // Check what has been sent
  if (!exist(rotx)) { rotx = 60; }
  if (!exist(rotz)) { rotz = 30; }
  if (!exist(scale)) { scale = 1; }
  if (!exist(scalez)) { scalez = 1; }

  if (!exist (N)) { N = DefPlotWin; }

  // Check the existence of p.[N]
  if (!exist(p.[N])) 
  { 
    pobj_Create (N); 
  }
    
  fprintf (p.[N].prog, "set view %f,%f,%f,%f\n", ...
           rotx, rotz, scale, scalez);
};

autotics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set xtics auto\n");
    fprintf (p.[N].prog, "set ytics auto\n");
    fprintf (p.[N].prog, "set ztics auto\n");
#     fprintf (p.[N].prog, "set nox2tics\n");
#     fprintf (p.[N].prog, "set noy2tics\n");
#     fprintf (p.[N].prog, "set nomxtics\n");
#     fprintf (p.[N].prog, "set nomytics\n");
#     fprintf (p.[N].prog, "set nomztics\n");
#     fprintf (p.[N].prog, "set nomx2tics\n");
#     fprintf (p.[N].prog, "set nomy2tics\n");
};

xtics = function ( start, incr, end, N )
{
    //Four arguments
    if (exist(N)) {
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set xtics %s,%s,%s\n",start,incr,end);
    }
    // No arguments
    if (!exist(start)) { 
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set xtics\n");
    }
    // One argument
    if (exist(start) && !exist(incr)) { 
	if (class(start) == "num") { // start is the plot window or a number
	    N = start;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set xtics\n");
	    else // start is a string
	    N = DefPlotWin;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set xtics %s\n",start);
	}
    }
    // Two arguments
    if (exist(incr) && !exist(end)) {
	if (class(start) == "string") { // start is a string
	    N = incr;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set xtics %s\n",start);
	    else
	    if (max(size(start)) > 1) { // start is a vector
		N = incr;
		if (!exist(p.[N])) { pobj_Create (N); }
		string = "set xtics (" + num2str(start[1])[1];
		for (i in 2:max(size(start))) {
		    if(class(i) == "string"){
			string = string + "," + i;
			else
			string = string + "," + num2str(start[i])[1];
		    }
		}
		else // start is not a vector
		N = DefPlotWin;
		if (!exist(p.[N])) { pobj_Create (N); }
		fprintf (p.[N].prog, "set xtics %s,%s\n",start,incr);
	    }
	}
    }
    // Three arguments
    if (exist(end) && !exist(N)) {
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set xtics %s,%s,%s\n",start,incr,end);
    }
};

noxtics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set noxtics\n");
};

x2tics = function ( start, incr, end, N )
{
    //Four arguments
    if (exist(N)) {
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set x2tics %s,%s,%s\n",start,incr,end);
    }
    // No arguments
    if (!exist(start)) { 
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set x2tics\n");
    }
    // One argument
    if (exist(start) && !exist(incr)) { 
	if (class(start) == "num") { // start is the plot window
	    N = start;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set x2tics\n");
	    else // start is a string
	    N = DefPlotWin;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set x2tics %s\n",start);
	}
    }
    // Two arguments
    if (exist(incr) && !exist(end)) {
	if (class(start) == "string") { // start is a string
	    N = incr;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set x2tics %s\n",start);
	    else
	    if (max(size(start)) > 1) { // start is a vector
		N = incr;
		if (!exist(p.[N])) { pobj_Create (N); }
		string = "set x2tics (" + num2str(start[1])[1];
		for (i in 2:max(size(start))) {
		    if(class(i) == "string"){
			string = string + "," + i;
			else
			string = string + "," + num2str(start[i])[1];
		    }
		}
		else // start is not a vector
		N = DefPlotWin;
		if (!exist(p.[N])) { pobj_Create (N); }
		fprintf (p.[N].prog, "set x2tics %s,%s\n",start,incr);
	    }
	}
    }
    // Three arguments
    if (exist(end) && !exist(N)) {
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set xytics %s,%s,%s\n",start,incr,end);
    }
};

nox2tics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nox2tics\n");
};

ytics = function ( start, incr, end, N )
{
    //Four arguments
    if (exist(N)) {
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set ytics %s,%s,%s\n",start,incr,end);
    }
    // No arguments
    if (!exist(start)) { 
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set ytics\n");
    }
    // One argument
    if (exist(start) && !exist(incr)) { 
	if (class(start) == "num") { // start is the plot window
	    N = start;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set ytics\n");
	    else // start is a string
	    N = DefPlotWin;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set ytics %s\n",start);
	}
    }
    // Two arguments
    if (exist(incr) && !exist(end)) {
	if (class(start) == "string") { // start is a string
	    N = incr;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set ytics %s\n",start);
	    else
	    if (max(size(start)) > 1) { // start is a vector
		N = incr;
		if (!exist(p.[N])) { pobj_Create (N); }
		string = "set ytics (" + num2str(start[1])[1];
		for (i in 2:max(size(start))) {
		    if(class(i) == "string"){
			string = string + "," + i;
			else
			string = string + "," + num2str(start[i])[1];
		    }
		}
		else // start is not a vector
		N = DefPlotWin;
		if (!exist(p.[N])) { pobj_Create (N); }
		fprintf (p.[N].prog, "set ytics %s,%s\n",start,incr);
	    }
	}
    }
    // Three arguments
    if (exist(end) && !exist(N)) {
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set ytics %s,%s,%s\n",start,incr,end);
    }
};

noytics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set noytics\n");
};

y2tics = function ( start, incr, end, N )
{
    //Four arguments
    if (exist(N)) {
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set y2tics %s,%s,%s\n",start,incr,end);
    }
    // No arguments
    if (!exist(start)) { 
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set y2tics\n");
    }
    // One argument
    if (exist(start) && !exist(incr)) { 
	if (class(start) == "num") { // start is the plot window
	    N = start;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set y2tics\n");
	    else // start is a string
	    N = DefPlotWin;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set y2tics %s\n",start);
	}
    }
    // Two arguments
    if (exist(incr) && !exist(end)) {
	if (class(start) == "string") { // start is a string
	    N = incr;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set y2tics %s\n",start);
	    else
	    if (max(size(start)) > 1) { // start is a vector
		N = incr;
		if (!exist(p.[N])) { pobj_Create (N); }
		string = "set y2tics (" + num2str(start[1])[1];
		for (i in 2:max(size(start))) {
		    if(class(i) == "string"){
			string = string + "," + i;
			else
			string = string + "," + num2str(start[i])[1];
		    }
		}
		else // start is not a vector
		N = DefPlotWin;
		if (!exist(p.[N])) { pobj_Create (N); }
		fprintf (p.[N].prog, "set y2tics %s,%s\n",start,incr);
	    }
	}
    }
    // Three arguments
    if (exist(end) && !exist(N)) {
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set y2tics %s,%s,%s\n",start,incr,end);
    }
};

noy2tics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set noy2tics\n");
};

ztics = function ( start, incr, end, N )
{
    if (!exist(start)) { start = 0; }
    if (!exist(p.[start])) { pobj_Create (start); }
    if (!exist(end)) {
	if (max(size(start)) > 1) {
	    if (!exist(incr)) { incr = 0; }
	    if (!exist(p.[incr])) { pobj_Create (incr); }
	    string = "set ztics (" + num2str(start[1])[1];
	    for (i in 2:max(size(start))) {
		string = string + "," + num2str(start[i])[1];
	    }
	    string = string + " nomirror)\n";
	    fprintf (p.[incr].prog, string);
	    else
	    if (!exist(incr)) {
		fprintf (start, p.[start].prog, "set ztics nomirror\n");
		else
		fprintf (start, p.[start].prog, "set ztics %s,%s nomirror\n",start,incr);
	    }
	}
	else
	// Set the default plot-object
	if (!exist(N)) { N = DefPlotWin; }
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set ztics %s,%s,%s nomirror\n",start,incr,end);
    }
};

noztics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set noztics\n");
};

mxtics = function ( freq, N )
{
    if (!exist(freq)) { error("No value given!"); }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set mxtics %s\n",freq);
};

nomxtics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nomxtics\n");
};

mx2tics = function ( freq, N )
{
    if (!exist(freq)) { error("No value given!"); }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set mx2tics %s\n",freq);
};

nomx2tics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nomx2tics\n");
};

mytics = function ( freq, N )
{
    if (!exist(freq)) { error("No value given!"); }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set mytics %s\n",freq);
};

nomytics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nomytics\n");
};

my2tics = function ( freq, N )
{
    if (!exist(freq)) { error("No value given!"); }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set my2tics %s\n",freq);
};

nomy2tics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nomy2tics\n");
};

mztics = function ( freq, N )
{
    if (!exist(freq)) { error("No value given!"); }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set mztics %s\n",freq);
};

nomztics = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nomztics\n");
};

xtr = function ( start, incr, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    xrange( start, end, N );
    xtics( start, incr, end, N );
};

ytr = function ( start, incr, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    yrange( start, end, N );
    ytics( start, incr, end, N );
};

ztr = function ( start, incr, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    zrange( start, end, N );
    ztics( start, incr, end, N );
};

x2tr = function ( start, incr, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    x2range( start, end, N );
    x2tics( start, incr, end, N );
};

y2tr = function ( start, incr, end, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    y2range( start, end, N );
    y2tics( start, incr, end, N );
};

loglog = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set logscale xy\n");
};

semilogx = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale xy\n");
    fprintf (p.[N].prog, "set logscale x\n");
};

semilogy = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale xy\n");
    fprintf (p.[N].prog, "set logscale y\n");
};

nolog = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale xy\n");
};

loglog2 = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set logscale x2\n");
    fprintf (p.[N].prog, "set logscale y2\n");
};

semilogx2 = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale y2\n");
    fprintf (p.[N].prog, "set logscale x2\n");
};

semilogy2 = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale x2\n");
    fprintf (p.[N].prog, "set logscale y2\n");
};

nolog2 = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nologscale x2\n");
    fprintf (p.[N].prog, "set nologscale y2\n");
};


autoscale = function ( axis, N )
{
    if (!exist(axis) && !exist(N)) {
	N = DefPlotWin;
	axis = ["xy","x2","y2"];
    }

    if (exist(axis) && !exist(N)) {
	if (class(axis) == "num") {
	    N = axis;
	    axis = ["xy","x2","y2"];
	}
    }

    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (class(axis) != "string") {
	error ("Must be of class string");
    }
    for (i in 1:max(size(axis))) {
	fprintf (p.[N].prog, "set autoscale %s\n",axis[i]);
    }
};

noautoscale = function ( axis, N )
{
    if (!exist(axis) && !exist(N)) {
	N = DefPlotWin;
	axis = ["xy","x2","y2"];
    }
    if (exist(axis) && !exist(N)) {
	if (class(axis) == "num") {
	    N = axis;
	    axis = ["xy","x2","y2"];
	}
    }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    if (class(axis) != "string") {
	error ("Must be of class string");
    }
    for (i in 1:max(size(axis))) {
	fprintf (p.[N].prog, "set noautoscale %s\n",axis[i]);
    }
};

psize = function ( x_size, y_size, N )
{
    // Three arguments
    if (exist(N)) { 
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set size %s,%s\n",x_size,y_size);
    }
    // No arguments
    if (!exist(x_size)) { 
	N = DefPlotWin;
	if (!exist(p.[N])) { pobj_Create (N); }
	fprintf (p.[N].prog, "set size\n");
    }
    // One argument
    if (exist(x_size) && !exist(y_size)) {
	if (class(x_size) == "string") {
	    N = DefPlotWin;
	    if (!exist(N)) { N = DefPlotWin; }
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set size %s\n",x_size);
	    else
	    N = x_size;
	    if (!exist(N)) { N = DefPlotWin; }
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set size\n");
	}
    }
    // Two arguments
    if (exist(y_size) && !exist(N_size)) {
	if (class(x_size) == "string") {
	    N = y_size;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set size %s\n",x_size);
	    else
	    N = DefPlotWin;
	    if (!exist(p.[N])) { pobj_Create (N); }
	    fprintf (p.[N].prog, "set size %s,%s\n",x_size,y_size);
	}
    }
};

origin = function ( x_origin, y_origin, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    
    fprintf (p.[N].prog, "set origin %s,%s\n",x_origin,y_origin);
};

xlabel = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set xlabel \"%s\"\n",string);
};

ylabel = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set ylabel \"%s\"\n",string);
};

zlabel = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set zlabel \"%s\"\n",string);
};

x2label = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set x2label \"%s\"\n",string);
};

y2label = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set y2label \"%s\"\n",string);
};

title = function ( string, N )
{
    // Just to be safe.
    if (!exist (string)) { string = ""; }

    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }

    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    fprintf (p.[N].prog, "set title \"%s\"\n",string);
};

//
// For back compatibility.
//

ptitle = title;


notitle = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set title\n");
};

linestyle = function ( string, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }

    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    if (!exist (string)) 
    {
      string = ["lines", "lines", "lines", "lines", "lines", "lines",...
                "lines", "lines", "lines", "lines", "lines", "lines"];
    }

    p.[N].linetype = string;
};

//
// Get the current linestyle/linetype...
//

get_linetype = function ( N,  I )
{
  //
  // Check the array for a specified linestyle.
  // If not found, then provide a default value
  // (lines).
  //

  if (length (length (p.[N].linetype)) >= I)
  {
    return p.[N].linetype[I];
  else
    return "lines";   # The default
  }
};


linenumbers = function ( nrs, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if(exist(nrs)){
	p.[N].lnr = nrs;
	p.[N].pnr = 0;
	else
	p.[N].lnr = [1,2,3,4,5,6,7,8];
	p.[N].pnr = 0;
    }
};

pointsize = function ( size, N )
{
    // Set the default plot-object.
    if (!exist(N)) { N = DefPlotWin; }

    // Default pointsize.
    if (!exist (size)) { size = 1; }

    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    // Set the pointsize.
    fprintf (p.[N].prog, "set pointsize %s\n",size);
};

pevery = function ( nr, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    p.[N].every = nr;
};

pformat = function ( string1, string2, N )
{
    // Check/set the default format
    if (!exist (string1)) { string1 = "%g"; }

    if (!exist(string2)) {
	string2 = "xy";
    }
    if (exist(string2) && !exist(N)) {
	if (class(string2) == "num") {
	    N = string2;
	    string2 = "xy";
	}
    }

    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }

    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    fprintf (p.[N].prog, "set format %s \"%s\"\n",string2,string1);
};

key = function ( x, y, z, N )
{
    if (exist(x) && !exist(y)) {
	if (class(x) == "num") {
	    N = x;
	    else
	    just = x;
	}
    }
    if (exist(x) && exist(y)) {
	if (class(x) == "string") {
	    just = x;
	    N = y;
	}
    }
    if (!exist(z)) { z = 0; }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (exist(just)) {
	fprintf (p.[N].prog, "set key %s\n",just);
	else
	if (exist(y)) {
	    fprintf (p.[N].prog, "set key right\n");
//	    fprintf (p.[N].prog, "set key %s,%s,%s\n",x,y,z);
	    else
	    fprintf (p.[N].prog, "set key right\n");
	}
    }
};

nokey = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nokey\n");
};

grid = function ( string, N )
{
    if (!exist(string)){
	if (class(string) == "num"){
	    N = string;
	    clear(string);
	}
    }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (exist(string)){
	fprintf (p.[N].prog, "set grid %s\n",string);
	else
	fprintf (p.[N].prog, "set grid\n");
    }
};

nogrid = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set nogrid\n");
};

zeroaxis = function ( axis, N )
{
    if (!exist(axis) && !exist(N)) {
	N = DefPlotWin;
	axis = "";
    }
    if (exist(axis) && !exist(N)) {
	if (class(axis) == "num") {
	    N = axis;
	    axis = "";
	}
    }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    if (class(axis) != "string") {
	error ("Must be of class string");
    }
    for (i in 1:max(size(axis))) {
	fprintf (p.[N].prog, "set %szeroaxis\n",axis[i]);
    }
};

nozeroaxis = function ( axis, N )
{
    if (!exist(axis) && !exist(N)) {
	N = DefPlotWin;
	axis = "";
    }
    if (exist(axis) && !exist(N)) {
	if (class(axis) == "num") {
	    N = axis;
	    axis = "";
	}
    }
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }

    if (class(axis) != "string") {
	error ("Must be of class string");
    }
    for (i in 1:max(size(axis))) {
	fprintf (p.[N].prog, "set no%szeroaxis\n",axis[i]);
    }
};

label = function ( label, x, y, z, labelnr, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (class(label) != "string") {
	error ("The label must be a string!")
    }
    if (exist(labelnr)) {
	fprintf (p.[N].prog, "set label %s \"%s\" at %s,%s norot\n", labelnr, label, x, y, z);
	else
	fprintf (p.[N].prog, "set label \"%s\" at %s,%s norot\n", label, x, y, z);
    }
};

nolabel = function ( labelnr, N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    if (exist(labelnr)) {
	fprintf (p.[N].prog, "set nolabel %i\n", labelnr);
	else
	fprintf (p.[N].prog, "set nolabel\n");
    }
};

ptime = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set time\n");
};

noptime = function ( N )
{
    // Set the default plot-object
    if (!exist(N)) { N = DefPlotWin; }
    // Check the existence of p.[N]
    if (!exist(p.[N])) { pobj_Create (N); }
    fprintf (p.[N].prog, "set notime\n");
};


//====================
// Static Functions	=
//====================

//
// Plot a list.
//

plotl = function ( data, N, keynames )
{
  I = 1;
  for (i in members (data)) 
  {
    if (class(data.[i]) == "num") 
    {
      // Check for plotting keynames
      plotm (data.[i], N, I, keynames);
      I++;
    else
      if (class(data.[i]) == "string") 
      {
	if (strsplt(data.[i])[1] == "w") 
	{
	  string = "";
	  for (j in 6:length(strsplt(data.[i]))) 
	  {
	    string = string + strsplt(data.[i])[j];
	  }
	  p.[N].linetype = string;
	else
	  if (strsplt(data.[i])[1] == "t") 
	  {
	    string = "";
	    for (j in 7:length(strsplt(data.[i]))) 
	    {
	      string = string + strsplt(data.[i])[j];
	    }
	    p.[N].keyname = string;
	  }
	}
      }
    }
  }
};



//
// Set-Up to plot a matrix. Columns 2...N against 1st column.
//

plotm = function ( data, N, I, keynames )
{
  // Determine how many lines to draw
  nplot = max([1, data.nc - 1]);

  if (nplot > data.nr) 
  {
    printf (" Plot %i columns, are you sure [y(es)/n(o)/t(ranspose)" + ...
            " (def. no)] ? ", data.nc);
    ans = getline ("stdin");
    if (ans.[1] != "y" && ans.[1] != "t") {
      return 0;
    }
    if (ans.[1] == "t") {
      data = data';
      nplot = max([1, data.nc - 1]);
    }
  }

  // Generate two column matrices for plot program
  p.[N].repnr = 1;
  for (i in 1:nplot) 
  {
    // Create tmp-file-name
    fn = pobj_TmpFileName (N, i, I);
    
    // Check for plotting keynames
    p.[N].lpnr = p.[N].lpnr + 1;
    pobj_PlotKeyName (keynames, p.[N].lpnr, N );
    
    // Write data to tmp-file
    // Add to plot command

    if (nplot == 1) 
    {
      pobj_WriteData (real(data), fn);
      pobj_PlotCmd (N, fn, i+I-1);
    else
      pobj_WriteData (real(data[;1,i+1]), fn);
      pobj_PlotCmd (N, fn, i+I-1);
    }
  }
};


replotl = function ( data, N, keynames )
{
    I = 1;
    for (i in members (data)) {
	if (class(data.[i]) == "num") {
	    // Check for plotting keynames
	    pobj_PlotKeyName (keynames, I, N);
	    replotm (data.[i], N, I);
	    I++;
	}
    }
};


replotm = function ( data, N, I, keynames )
{
    // Determine how many lines to draw
    nplot = max([1, data.nc - 1]);

    if (nplot > data.nr) {
	printf (" Plot %i columns, are you sure [y(es)/n(o)/t(ranspose) (def. no)] ? ", data.nc);
	ans = getline ("stdin");
	if (ans.[1] != "y" && ans.[1] != "t") {
	    return 0;
	}
	if (ans.[1] == "t") {
	    data = data';
	    nplot = max([1, data.nc - 1]);
	}
    }

    // Generate two column matrices for plot program
    for (i in 1:nplot) {
	// Create tmp-file-name
	fn = pobj_TmpFileName (N, i, I);

	// Check for plotting keynames
	p.[N].lpnr = p.[N].lpnr + 1;
	pobj_PlotKeyName (keynames, p.[N].lpnr, N);

	// Write data to tmp-file
	// Add to plot command
	if (nplot == 1) {
	    pobj_WriteData (real(data), fn);
	    pobj_RePlotCmd (N, fn, i+I-1);
	    else
	    pobj_WriteData (real(data[;1,i+1]), fn);
	    pobj_RePlotCmd (N, fn, i+I-1);
	}
    }
};



//
// Set-Up to splot a matrix. 
//

splotm = function ( datax, datay, datazz, nsplot, N, I, keynames )
{
    // Generate three column matrices for splot program
    p.[N].repnr = 1;
    // Write data to new data-file
    // Create tmp-file-name
    for (i in 1:nsplot) {
	fn = pobj_TmpFileName (N, i, I);

	// Check for plotting keynames
	pobj_PlotKeyName (keynames, I, N);
	newdata = [];
	for (m in 1:datax.nr){
	    rownumber = 0;
	    for (n in 1:datay.nr){
		rownumber = rownumber+1;
		// x-value
		newdata[rownumber;1] = datax[m];
		// y-value
		newdata[rownumber;2] = datay[n];
		// z-value
		newdata[rownumber;3] = datazz[m;n+(i-1)*datax.nr];
# 		newdata[rownumber;3] = datazz[n+(i-1)*datax.nr;m];
	    }
# 	    writem (fn,real(newdata[;1,2,3]));
# 	    writem (fn,"");
#	    // Add an empty line
# 	    close (fn);
# 	    system(" echo '' >> " + fn);
# 	    open (fn,"append");
	    a=real(newdata[;1,2,3]);
	    for(row in 1:size(a)[1]) {
		for(column in 1:size(a)[2]) {
		    fprintf(fn, a[row;column]);
		    fprintf(fn, "\t");
		}
		fprintf(fn, "\n");
	    }
	    fprintf(fn, "\n");
	}
	close (fn);
	// Add to splot command
	pobj_SPlotCmd (N, fn, i+I-1);
    }
};

splotm_cyl = function ( datara, datath, datazz, nsplot, N, I, keynames )
{
    // Generate three column matrices for splot program
    p.[N].repnr = 1;
    // Write data to new data-file
    // Create tmp-file-name
    for (i in 1:nsplot) {
	fn = pobj_TmpFileName (N, i, I);

	// Check for plotting keynames
	pobj_PlotKeyName (keynames, I, N);
	newdata = [];
	for (m in 1:datath.nr){
	    rownumber = 0;
	    for (n in 1:datara.nr){
		rownumber = rownumber+1;
		// Theta
		newdata[rownumber;1] = datath[m];
		// Radius
		newdata[rownumber;3] = datara[n];
		// z-value
		newdata[rownumber;2] = datazz[m;n+(i-1)*datara.nr];
# 		newdata[rownumber;2] = datazz[n+(i-1)*datara.nr;m];
	    }
# 	    writem (fn,real(newdata[;1,2,3]));
# 	    writem (fn,"");
#	    // Add an empty line
# 	    close (fn);
# 	    system(" echo '' >> " + fn);
# 	    open (fn,"append");
	    a=real(newdata[;1,2,3]);
	    for(row in 1:size(a)[1]) {
		for(column in 1:size(a)[2]) {
		    fprintf(fn, a[row;column]);
		    fprintf(fn, "\t");
		}
		fprintf(fn, "\n");
	    }
	    fprintf(fn, "\n");
	}
	close (fn);
	// Add to splot command
	pobj_SPlotCmd (N, fn, i+I-1);
    }
};


//
// Form a plain string to send to GNUPLOT as 
// command
//

plots = function ( data, N )
{
    // Send the string to GNUPLOT
    p.[N].pcmd = data + "\n";
};



//
// Create a plot-object.
//

pobj_Create = function ( N )
{
    // The tmp files to plot
    plist.files = "";
    // Where the plot command will go
    plist.pcmd = "";
    
    // Which y-axis
    plist.axes = FoS;
    // Which linestyle
    plist.linetype = ["lines", "lines", "lines", "lines", "lines", "lines",...
                      "lines", "lines", "lines", "lines", "lines", "lines"];
    // Every ? point
    plist.every = 1;
    // Which keykeyname
    plist.keyname = "";
    #for (i in 1:20) { plist.keyname[i] = "C-" + num2str(i); }
    // Number of plotted curves
    plist.pnr = 0;
    plist.lpnr = 0;
    plist.lnr = [1,2,3,4,5,6,7,8];
    // Number of replots
    plist.repnr = 1;
    // Some parameters to control the multiplot functions
    plist.mpp = [];       // Row, column, Tracking index, Plot Window N
    plist.mpo = [];       // Size and origin information
    
    // Init string for plotting program
    plist.init = "set nogrid\nset mxtics 5\nset mytics 5" + ...
                 "\nset data style lines\nset nokey" + ...
                 "\nset noxzeroaxis\nset noyzeroaxis\nset encoding iso_8859_1\n";

    // The program that draws the plot(s)
    sprintf (plist.prog, "|gnuplot -title 'RLaB2 #%s' #%i", N, N);

    // To remove tmp-files
#     sprintf (plist.rmf, "|rm -f `cat`");
    sprintf (plist.rmf, "|rm -f `cat` #%i", N);
    
    // Keep track of terminal type, and output
    plist.term = "X11";
    plist.output = "stdout";

    // Copy the local list into the static plot-object collection
    p.[N] = plist;
    fprintf (p.[N].prog, "%s", p.[N].init);
    fprintf (p.[N].prog, "set term %s\n", p.[N].term);
};



//
// Reset a plot object to plot new data.
// Close any existing tmp-files that belong
// to plot-object N. Reset the file-name list.
//

pobj_Reset = function ( N )
{
    pobj_Rm (N);
    p.[N].files = "";
    p.[N].pcmd = "";
    p.[N].pnr = 0;
};


//
// Destroy a plot-object
//

pobj_Destroy = function ( N )
{
    if (exist(p.[N])) {
	clear (p.[N]);	
    }
};



//
// Create a tmp-file name
//

pobj_TmpFileName = function ( N, i, j )
{
    // N = plot window
    // i = curve number
    // p.[N].repnr = number of replots
    // j = list number

    tmp = tmpnam();
    p.[N].files = p.[N].files + " " + tmp;
    return tmp;
};



//
// Add data-file (tmp-file) to plot object list of stuff
// to be plotted.
//

pobj_WriteData = function ( m, file )
{
  // Save and set the ASCII write format.
  sf = format (20,10);

  // Write the data.
  writem (file, m);

  // Restore the ASCII write format.
  format (sf);

  close (file);
};



//
// Create the command(s) to plot all the data in the plot object.
//

pobj_PlotCmd = function ( N, fn, I )
{
  if(p.[N].pnr >= length(p.[N].lnr))
  {
    p.[N].pnr = 1;
    else
    p.[N].pnr = p.[N].pnr + 1;
  }

  if(p.[N].pcmd == "") 
  {
    sprintf(tmp, "plot '%s' every %s axes %s title '%s' with %s %s", ...
            fn, p.[N].every, p.[N].axes, p.[N].keyname, ...
            get_linetype (N, I), p.[N].lnr[p.[N].pnr]);
  else
    sprintf(tmp, ", '%s' every %s axes %s title '%s' with %s %s", ...
            fn, p.[N].every, p.[N].axes, p.[N].keyname, ...
            get_linetype (N, I), p.[N].lnr[p.[N].pnr]);
  }

  p.[N].pcmd = p.[N].pcmd + tmp;
  p.[N].keyname = "";

};



//
// Create the command(s) to replot all the data in the plot object.
//

pobj_RePlotCmd = function ( N, fn, I )
{
    if(p.[N].pnr >= length(p.[N].lnr)){
	p.[N].pnr = 1;
	else
	p.[N].pnr = p.[N].pnr + 1;
    }
    if(p.[N].pcmd == "") {
	sprintf(tmp, "replot '%s' every %s axes %s title '%s' with %s %s", ...
                fn, p.[N].every, p.[N].axes, p.[N].keyname, ...
                get_linetype (N, I), p.[N].lnr[p.[N].pnr]);
	else
	sprintf(tmp, ", '%s' every %s axes %s title '%s' with %s %s", ...
                fn, p.[N].every, p.[N].axes, p.[N].keyname, ...
                get_linetype (N, I), p.[N].lnr[p.[N].pnr]);
    }
    p.[N].pcmd = p.[N].pcmd + tmp;
    p.[N].keyname = "";
};


//
// Create the command(s) to splot all the data in the splot object.
//

pobj_SPlotCmd = function ( N, fn, I )
{
    if(p.[N].pnr >= length(p.[N].lnr)){
	p.[N].pnr = 1;
	else
	p.[N].pnr = p.[N].pnr + 1;
    }
    if(p.[N].pcmd == "") {
	sprintf(tmp, "splot '%s' every %s title '%s' with %s %s", ...
                fn, p.[N].every, p.[N].keyname, get_linetype (N, I), ...
                p.[N].lnr[p.[N].pnr]);
	else
	sprintf(tmp, ", '%s' every %s title '%s' with %s %s", ...
                fn, p.[N].every, p.[N].keyname, get_linetype (N, I), ...
                p.[N].lnr[p.[N].pnr]);
    }

    p.[N].pcmd = p.[N].pcmd + tmp;
    p.[N].keyname = "";
};



//
// Check for plotting keynames
//

pobj_PlotKeyName = function (keynames, i, N)
{
    if (exist(keynames)) {
	if (i > max(size(keynames))) {
	    p.[N].keyname = "";
	    else
	    p.[N].keyname = keynames[i];
	}
    }
#    if (!exist(keynames)) {
#	sprintf (p.[N].keyname, "%.4g", i+1);
#	p.[N].keyname = "C-" + p.[N].keyname;
#	else
#	if (i > max(size(keynames))) {
#	    sprintf (p.[N].keyname, "%.4g", i+1);
#	    p.[N].keyname = "C-" + p.[N].keyname;
#	    else
#	    p.[N].keyname = keyname!s[i];
#	}
#    }
};



//
// Force the plot program to create the plot
// Other functions set up gnuplot, and the temporary
// data files. This function produces the final plot 
// command.
//

pobj_Plot = function ( N, debug )
{
  if (!exist (debug)) { debug = 0; }

  if (debug)
  {
    fprintf ("stdout", "%s\n", p.[N].pcmd);
  }
  
  fprintf (p.[N].prog, "%s\n\n", p.[N].pcmd);
};



//
// Setup so the tmp-files can be removed
// The plot-object element rmf contains a command
// that looks like: rm -f `cat`
//

pobj_SetRm = function ( N, debug )
{
  if (!exist (debug)) { debug = 0; }

  if (debug)
  {
    fprintf("stdout", "%s", p.[N].files);
  else
    fprintf(p.[N].rmf, "%s", p.[N].files);
  }
};



//
// Remove the tmp-files
// Closing the process forces execution of the
// rm command (on Unix systems).
//

pobj_Rm = function ( N )
{
  if (length(p.[N].files) != 0) 
  {
    close (p.[N].rmf);
  }
};
