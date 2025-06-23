savepdf2 = function(E, SN ){
# savepdf2 = function(plt1, plt2, SN ){
// based on gnuplot_eps.diary.r
//
// saves pdf plot copy to hardcoded /var/www/html directory
//   E -> list containnig four plots
//   plt1 -> list containnig first plot
//   plt2 -> list containnig second plot
//   SN   -> string containing serial number (to identify saved plot)
//   
//
//
///////////////////////////////////////////////////
# require libgnuplot
plt1 = E.plt2;
plt2 = E.plt3;
plt3 = E.plt4;
plt4 = E.plt5;

# gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot5junk.pdf'");           
gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot4junk.eps'");           
gnucmd("set multiplot");
gnucmd("set origin 0,0.55");
gnucmd("set size 1, 0.45");
# gnuxlabel("Time [s]");
gnuxlabel("Travel Distance [mm]");
#gnuylabel("Travel [mm]");
gnuylabel("Error [um]");
# gnulegend(["Measured travel","Ideal"]);
gnulegend("(First motion) FWD Travel Error [um]");
#gnutitle("Coarse motion vs. constant slope "+SN);
gnutitle("Deviation from linear motion of FH-MAA "+SN);
gnuformat ("with lines");
gnulimits(-1.5,9,,,,);
gnuplot(plt1);
 // plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size 1, 0.45");
# gnuxlabel("Time [s]");
gnuxlabel("Travel Distance [mm]");
#gnulegend("Travel Error [um]");
gnulegend("(Second motion) REV Travel Error [um]");
gnutitle("Deviation from linear motion of FH-MAA "+SN);
gnuylabel("Error [um]");
gnulimits(-1.5,9,,,,);
gnuplot(plt2);
gnuclose(1);
# if(dir == +1) {UD = "_FWD";}
# if(dir == -1) {UD = "_REV";}
UD="";
fn2=SN + UD + ".pdf";
# RPi version
# system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/var/www/html/"+fn2);
system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/tmp/"+fn2);
######################################################################

# gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot5junk.pdf'");           
gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot34junk.eps'");           
gnucmd("set multiplot");
gnucmd("set origin 0,0.55");
gnucmd("set size 1, 0.45");
# gnuxlabel("Time [s]");
gnuxlabel("Travel Distance [mm]");
#gnuylabel("Travel [mm]");
gnuylabel("Error [um]");
# gnulegend(["Measured travel","Ideal"]);
gnulegend("(Third motion) FWD Travel Error [um]");
#gnutitle("Coarse motion vs. constant slope "+SN);
gnutitle("Deviation from linear motion of FH-MAA "+SN);
gnuformat ("with lines");
gnulimits(-1.5,9,,,,);
gnuplot(plt3);
 // plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size 1, 0.45");
# gnuxlabel("Time [s]");
gnuxlabel("Travel Distance [mm]");
#gnulegend("Travel Error [um]");
gnulegend("(Fourth motion) REV Travel Error [um]");
gnutitle("Deviation from linear motion of FH-MAA "+SN);
gnuylabel("Error [um]");
gnulimits(-1.5,9,,,,);
gnuplot(plt4);
gnuclose(1);
# if(dir == +1) {UD = "_FWD";}
# if(dir == -1) {UD = "_REV";}
UD="";
fn34=SN + UD + "34.pdf";
# RPi version
# system("cat /tmp/gnuplot34junk.eps |/usr/bin/epstopdf -f -o=/var/www/html/"+fn34);
system("cat /tmp/gnuplot34junk.eps |/usr/bin/epstopdf -f -o=/tmp/"+fn34);
######################################################################
sleep(2);
system("/usr/bin/pdfunite /tmp/"+fn2+" /tmp/"+fn34+" /var/www/html/"+fn2);

};
