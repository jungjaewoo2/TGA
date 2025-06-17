//
//
//
if (isfile("plplot.r"))
{
  rfile plplot
}
if (length(plwins().available)<1)
{
  plwins (1,"xwin",[1200,960],[3000,0]);
}

x  = [0:4:1/32]';
dx = 0.1 * gaussian(x);

z1 = gaussian(1000,1);
z1_hist = hist(z1,x);

z2 = uniform(1000,1);
z2_hist = hist(z2,x);

y1 = pi * sin(pi * x);
dy1 = 0.1 * gaussian(x) .* y1;

_pl_data = << ...
  a1 = z1_hist;
  a2 = z2_hist;
  a3=[x,y1]; ...
  a4=[x,y1,x-dx,x+dx,y1-dy1,y1+dy1]>>;

tic();
i=0;
while(toc()<1)
{
  plwin  (1);
  pltitle ("looking for strings",0.75);
  plimits(0,4,-4,4);
  plxlabel("X-coordinate");
  plylabel("Y-coordinate");
  plegend(["histogram 1", "histogram 2", "uh-oh", "with errors on top"], 0.5,"rits",[0.1,0.15]);
  plformat([ ...
      "with lines lt 1 lc rgb black",...
      "with lines lt 1 lc rgb red",...
          "with lines lt 1 lc rgb brown", ...
              "with xyerrorbars lt 1 pc rgb green pt 0 ps 1", ...
                  []]);
  plxtics (0.5,5);
  plytics (1,5);
  plplot  (_pl_data);
  i++;
  sleep(0.002);
}
_plprint ("eg6.eps", "psc");
printf("PLPLOT: In 1 second 'plplot()' was executed %g times\n", i);



