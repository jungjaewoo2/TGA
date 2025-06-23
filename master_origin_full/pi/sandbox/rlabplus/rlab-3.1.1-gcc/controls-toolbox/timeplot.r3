//-----------------------------------------------------------------------
//
// timeplot
//
// Syntax: timeplot(iy,iu,f,ntim,Y);
//
// This routine plot the time response data in the block data
// matrix Y. The inputs are:
//
//    iy, iu -- Output number and input number to plot
//    f -- Data sample rate
//    ntim -- Number of time blockes stored in Y
//    Y -- Block-storage data matrices
//    name1,name2 -- Annotation labels for the plots
//
//
// Originally written by Lee D. Peterson for MATLAB
// Modified and ported to RLaB be Jeffrey B. Layton
// Version JBL 940519
//-----------------------------------------------------------------------

timeplot = function(iy,iu,f,ntim,Y,name)
{
   global (_rlab_config,_ctb2_window)

   if (nargs != 6) {
       error("timeplot: Wrong number of input arguments");
   }


// Get the sizes of the blocks in Y
   m=Y.nc;
   l=Y.nr/ntim;
   if (fix(l) != l) {
       error("Y dimension is not an even multiplier of given ntim.");
   }

// Check the input / output designators
   printf("iy = %i \n",iy);
   printf("iu = %i \n",iu);
   printf("name = %s \n",name);
   if ( (iy < 0) || (iy > l) || (iu < 0) || (iu > m) ) {
       error("Selected plot is out-of-bounds.");
   }


   if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
   {
     if (exist(_ctb2_window))
     {
        plwin(_ctb2_window);
     } else {
        _ctb2_window = plstart();
     }
   }
// Set-up the time axis vector
   t=[0:(ntim-1)/f:1/f]';

// Label the axes
   xlabel("t (seconds)");

// Title the plot
   strout="Response y"+int2str(iy)+" to Input u"+int2str(iu);
   printf("%s",strout);
   printf(" \n");
   pltitle(strout);

// Label the data
   plegend(name);

// Make the plot
   yindex=iy:l*ntim:l;

   plot([t,Y[yindex;iu]]);
   xlabel("");
   
};

