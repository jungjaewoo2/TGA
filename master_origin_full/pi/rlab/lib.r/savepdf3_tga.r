savepdf3_tga = function(E, SN ){
	#HACKED FOR TESTING DO NOT USE
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

# replace spaces with underscores
SN=gsub("_"," ",SN).string;

# gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot5junk.pdf'");           
gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot4junk.eps'");           
gnucmd("set multiplot");
gnucmd("set origin 0,0.55");
gnucmd("set size 1, 0.45");
gnuxlabel("Travel Distance [mm]");
gnuylabel("Error [um]");
gnulegend("(First motion) FWD Travel Error [um]");
gnutitle("Deviation from linear motion of HH-MAA "+SN);
gnuformat ("with lines");
gnulimits(-0.5,8.5,,,,);
gnuplot(plt1);
 // plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size 1, 0.45");
gnuxlabel("Travel Distance [mm]");
gnulegend("(Second motion) REV Travel Error [um]");
gnutitle("Deviation from linear motion of HH-MAA "+SN);
gnuylabel("Error [um]");
gnulimits(-8.5,0.5,,,,);
gnuplot(plt2);
gnuclose(1);
UD="";
fn2=SN + UD + ".pdf";
# RPi version
# system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/var/www/html/"+fn2);
system("cat /tmp/gnuplot4junk.eps |/usr/bin/epstopdf -f -o=/tmp/"+fn2);
######################################################################

gnuwins(1, "set term postscript enh color; set output '/tmp/gnuplot34junk.eps'");           
gnucmd("set multiplot");
gnucmd("set origin 0,0.55");
gnucmd("set size 1, 0.45");
gnuxlabel("Travel Distance [mm]");
gnuylabel("Error [um]");
gnulegend("(Third motion) FWD Travel Error [um]");
gnutitle("Deviation from linear motion of HH-MAA "+SN);
gnuformat ("with lines");
gnulimits(-0.5,8.5,,,,);
gnuplot(plt3);
 // plot no. 2 in the same device
gnucmd("set origin 0,0");
gnucmd("set size 1, 0.45");
gnuxlabel("Travel Distance [mm]");
gnulegend("(Fourth motion) REV Travel Error [um]");
gnutitle("Deviation from linear motion of HH-MAA "+SN);
gnuylabel("Error [um]");
gnulimits(-8.5,0.5,,,,);
gnuplot(plt4);
gnuclose(1);
UD="";
fn34=SN + UD + "34.pdf";
# RPi version
# system("cat /tmp/gnuplot34junk.eps |/usr/bin/epstopdf -f -o=/var/www/html/"+fn34);
system("cat /tmp/gnuplot34junk.eps |/usr/bin/epstopdf -f -o=/tmp/"+fn34);
######################################################################
sleep(2);
# system("/usr/bin/pdfunite /tmp/"+fn2+" /tmp/"+fn34+" /tmp/html/"+fn2);
system("/usr/bin/pdfunite /tmp/"+fn2+" /tmp/"+fn34+" /var/www/html/"+fn2);

};
