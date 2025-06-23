ap5 = function(Sig){
// Plot family of curves using autoplot
// Arguments:
//      Sig          - matrix to be plotted
// Optional Arguments:
//      dt           - increment (used only if column vector)
//      PlotTile     - specify or else assigns consecutive number
//      PlotFileName - specify or else creates temporary file in /var/tmp
//      PrintFlag    - "y" or "n"
//
// ToDo: fix printing

rfile tmpnam
global(WinNum);
local(pname,plotfile, apcmd, plt, x, M);


pname = "|ap -w ap" ;
plotfile = "/tmp/currentplot.ap" ;

plt = Sig;

#fprintf (plotfile, "%s\n", apcmd);
fprintf (plotfile,  "title Frequency response\n");
fprintf (plotfile,  "logx\n");
fprintf (plotfile,  "xscale Frequency [Hz]\n");
fprintf (plotfile,  "yscale Magnitude [dB]\n");
sleep 0.1;
write_ascii(plotfile,plt);
system("ap -w Magnitude " + plotfile + ";");
sleep 0.1;
#system("pd -w Magnitude " + plotfile + ";");
close(plotfile);

};
