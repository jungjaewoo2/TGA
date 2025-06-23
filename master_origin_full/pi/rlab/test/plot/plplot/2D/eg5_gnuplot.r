//
//
//

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
  gnuwins (1);
  gnuwin  (1);
  gnulimits(0,4,-4,4);
  gnuxlabel("X-coordinate");
  gnuylabel("Y-coordinate");
  gnulegend(["histogram 1", "histogram 2", "", "with errors on top"]);
  gnuformat([ ...
    "with lines lt 1 lc 16",...
    "with lines lt 1 lc 2",...
    "with lines lt 1 lc 1", ...
    "with xyerrorbars lt 1 pc 4 pt 1 ps 1", ...
  []]);
  gnuxtics (0.5,5);
  gnuytics (1,5);
  gnuplot  (_pl_data);
  i++;
}
printf("GNUPLOT: In 1 second 'gnuplot()' was executed %g times\n", i);



