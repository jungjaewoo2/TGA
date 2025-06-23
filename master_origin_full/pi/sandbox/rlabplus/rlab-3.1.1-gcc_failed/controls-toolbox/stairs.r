//------------------------------------------------------------------------
//
// stairs
//
// Syntax: A=stairs(x,y)
//
// This routine draws a stairstep graph of the elements in vector y
// at the locations specified in vector x. A stiarstep graph is similar
// to a bar graph, but the vertical lines dropping to the x-axis are
// omitted. Stairstep plots are useful for drawing time history plots
// of digital sampled-data systems.
//
// The routine can also be called as stairs(y), which draws a stairstep
// graph of the elements of vector y.
//
// Note: If the routine is called as stairs(x,y), then the values in x
//       must be evenly spaced in ascending order.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940906
//------------------------------------------------------------------------

stairs = function(x,y)
{
   global(_rlab_config,_ctb2_window)

   // Check dimensions
   if (nargs == 1) {
       n=length(x);
   else
       n=length(x);
       if (length(y) != n) {
           error("STAIRS:  Length(y) != Length(x)");
       }
   }
   ip=(n*2)-1;
   xp=zeros(ip,1);
   yp=zeros(ip,1);

   // Convert vectors x and y into arrays for plotting
   if (nargs == 2) {
       for (i in 1:n) {
            ip=(i-1)*2 + 1;
            xp[ip]=x[i];
            if (i != n) {
                xp[ip+1]=x[i+1];
            }
            yp[ip]=y[i];
            if (i != n) {
                yp[ip+1]=y[i];
            }
       }
   else
       for (i in 1:n) {
            ip=(i-1)*2 + 1;
            xp[ip]=i;
            if (i != n) {
                xp[ip+1]=i+1;
            }
            yp[ip]=x[i];
            if (i != n) {
                yp[ip+1]=x[i];
            }
       }
   }

   // Make the plot
   if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
   {
     if (exist(_ctb2_window))
     {
        plwin(_ctb2_window);
     else
        _ctb2_window = plstart();
     }
   }
   plot([xp,yp]);


   return << xp=xp; yp=yp >>;
};

