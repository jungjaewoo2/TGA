savepdf = function(plt1, plt2, SN, dir){
// based on gnuplot_eps.diary.r
//
// saves pdf plot copy to hardcoded /var/www/html directory
//   plt1 -> list containnig first plot
//   plt2 -> list containnig second plot
//   SN   -> string containing serial number (to identify saved plot)
//   dir  -> +1/-1 indicating direction of motion "forward" or "reverse" used to build a file name
//
//
///////////////////////////////////////////////////
# require libgnuplot

# gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot5junk.pdf'");           
gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot4junk.eps'");           
gnucmd("set multiplot");
gnucmd("set origin 0,0.55");
gnucmd("set size 1, 0.45");
gnuxlabel("Time [s]");
gnuylabel("Travel [mm]");
gnulegend(["Measured travel","Ideal"]);
gnutitle("Coarse motion vs. constant slope "+SN);
gnuformat ("with lines");
gnuplot(plt1);
 // plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size 1, 0.45");
gnuxlabel("Time [s]");
gnulegend("Travel Error [um]");
gnutitle("Deviation from linear motion of FH-MAA "+SN);
gnuylabel("Error [um]");
gnuplot(plt2);
gnuclose(1);
if(dir == +1) {UD = "_FWD";}
if(dir == -1) {UD = "_REV";}
fn2=SN + UD + ".pdf";
# RPi version
system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/var/www/html/"+fn2);
# system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/tmp/"+fn2);
};
