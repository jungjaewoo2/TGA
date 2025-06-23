ap2 = function(Sig,dt,PlotTile,PlotFileName,PrintFlag){
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

global(WinNum);
local(pname, plotfile, apcmd, plt, x, M);

if( exist(WinNum)) {
                     WinNum = WinNum + 1;
                  else
                     WinNum = 1;
}

if( exist(PlotTile)) {
                apcmd = "title "+PlotTile+"\n";
        else
                apcmd = "title ap"+num2str(WinNum)+"\n";

pname = "|ap -w ap" + num2str(WinNum);

if( (Sig.nr == 1) || (Sig.nc == 1)){
    // create x column
    if( exist(dt)) {
    x = (1:max(size(Sig)))*dt;
    else
    x = (1:max(size(Sig)));
    }
    x = x[:];
    Sig = Sig[:];

        if( max(imag(Sig)) == 0){
             plt=[x,Sig];
         else
             plt = [[x,real(Sig)];[x,imag(Sig)]];
             }

else
//    plt = Sig;
plt=[];
//M = Sig.nc;
for ( i in (2:Sig.nc)) { 
//x =Sig[;1,i];
plt= [plt;Sig[;1,i]];
}
}
if(exist(PlotFileName)) {
                 plotfile = PlotFileName;
                 else
                 plotfile = tmpnam() + ".ap";
}

fprintf (plotfile, "%s\n", apcmd);
write_ascii(plotfile,plt);
close(plotfile);

fprintf (pname, "%s\n", apcmd);
write_ascii(pname,plt);
close(pname);

#if(exist(PrintFlag)) {
#      if( PrintFlag == "y" ) {system("djet -a "+plotfile);
#      }}
if(!exist(PrintFlag)) {
// Should we print?
xstr.[1]="";
#printf( "Print the results? (y/n)" ); xstr = getline( "stdin" );
#if ( xstr.[1] == "y" )
#{
#system("djet -a "+plotfile);
#}
}
}
pause("Press return");
}
