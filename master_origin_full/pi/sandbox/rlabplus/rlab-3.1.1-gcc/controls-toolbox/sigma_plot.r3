sigma_plot = function (R)
{
   // plot Singular value frequency response of continuous linear systems.
   //
   //   call sequence:
   //
   //   R = sigma(a,b,c,d,w,invflag);
   //   sigma_plot(R);
   //
   global (_rlab_config,_ctb2_window)

   if (_rlab_config.plot_support=="pgplot"||_rlab_config.plot_support=="plplot")
   {
     if (exist(_ctb2_window))
     {
        plwin(_ctb2_window);
     } else {
        _ctb2_window = plstart();
     }
   }

   pltitle("Singular Value Frequency Response");
   xlabel("Frequency (rad/sec)");
   ylabel("Singular Values dB");
   plegend();
   plaxis("log");
   plot([R.w,20*log10(R.svout)]); 
   plaxis();
   plegend("default");
};

