//
// libhistogram.r
//
// by marijan kostrun, (c) 2000-2008

// contains:
//  ishist, ishist2
//  hist_line, hist2_surf

// are histograms normalized?
static (HISTOGRAM_NORM);
// HISTOGRAM_NORM = 0; // not: magnitude is the size of the bin
HISTOGRAM_NORM = 1; // yes: reduces to a probability distribution

//
// ishist: check whether a list is a histogram
//
ishist = function(x)
{
  if (class(x) != "list")
  { return 0; }

  if (all(members(x)!="bin"))
  { return 0; }

  if (all(members(x)!="range"))
  { return 0; }

  if (length(x.range)-length(x.bin)!=1)
  { return 0; }

  iok=1;
  i = 2;
  while(iok==1 && i<=length(x.range))
  {
    if(x.range[i] <= x.range[i-1])
    {
      iok = 0;
      break;
    }
    i++;
  }
  return iok;
};

//
// ishist2: check whether a list is a 2D-histogram
// <<bin; xrange; yrange>>
//
ishist2 = function(x)
{
  if (class(x) != "list")
  { return 0; }

  if (all(members(x)!="bin"))
  { return 0; }
  if (all(members(x)!="xrange"))
  { return 0; }
  if (all(members(x)!="yrange"))
  { return 0; }

  if (length(x.xrange)-(x.bin.nr)!=1)
  { return 0; }
  if (length(x.yrange)-(x.bin.nc)!=1)
  { return 0; }

  iok=1;
  i = 2;
  while(iok==1 && i<=length(x.xrange))
  {
    if(x.xrange[i] <= x.xrange[i-1])
    {
      iok = 0;
      break;
    }
    i++;
  }
  while(iok==1 && i<=length(x.yrange))
  {
    if(x.yrange[i] <= x.yrange[i-1])
    {
      iok = 0;
      break;
    }
    i++;
  }
  return iok;
};

// utility functions that determines plotting of the histograms
hist_norm = function ( i )
{
  if (exist(i))
  {
    if (i==0)
    {
      HISTOGRAM_NORM = 0;
    else
      HISTOGRAM_NORM = 1;
    }
  }

  return HISTOGRAM_NORM;
};

//
// utility functions for plotting libraries:
//  converts histogram data to a [x,y] line for
//  plotting in a standard graph window according to mode:
//    mode=1: bar plot, the width of a bar is the particular range
//    mode=2: line plot with the abscisa at the center of the range
hist_line = function (data, mode)
{
  if (!exist(mode))
  { mode = 1; }

  // data = <<bin;range>>
  if(!ishist(data))
  { return []; }

  s = 1;
  dn = sum(data.bin);
  if (mode==1)
  {
    //
    // mode=1: bar plot, the width of a bar is the particular range
    //
    npt = 3 * length(data.bin) + 1;
    p = zeros(npt, 2);
    p[1;1] = data.range[1];
    p[1;2] = 0;
    for (i in 1:length(data.bin))
    {
      // normalize histogram?
      if (HISTOGRAM_NORM)
      {  s = (data.range[i+1] - data.range[i]) * dn; }
      p[3*(i-1) + 2;1] = data.range[i];
      p[3*(i-1) + 2;2] = data.bin[i]/s;
      p[3*(i-1) + 3;1] = data.range[i+1];
      p[3*(i-1) + 3;2] = data.bin[i]/s;
      p[3*(i-1) + 4;1] = data.range[i+1];
      p[3*(i-1) + 4;2] = 0;
    }

  else
    //
    // mode=2: line plot with the abscisa at the center of the range
    //
    npt = length(data.bin);
    p   = zeros(npt, 2);
    for (i in 1:length(data.bin))
    {
      if (HISTOGRAM_NORM)
      {  s = (data.range[i+1] - data.range[i]) * dn; }
      p[i;1] = 0.5*(data.range[i]+data.range[i+1]);
      p[i;2] = data.bin[i]/s;
    }
  }

  return p;
};

hist2_surf = function (data)
{
  // data = <<bin;xrange;yrange>>
  if(!ishist2(data))
  { return []; }

  nx = data.bin.nr;
  ny = data.bin.nc;

  // construct centroid
  x = 0.5*( data.xrange[1:nx] + data.xrange[2:nx+1]);
  y = 0.5*( data.yrange[1:ny] + data.yrange[2:ny+1]);

  dx = data.xrange[2:nx+1] - data.xrange[1:nx];
  resize(dx,length(dx),1);    // always a column vector
  dy = data.yrange[2:ny+1] - data.yrange[1:ny];
  resize(dy,1,length(dy));    // always a row vector

  q = <<>>;
  q.x = x;
  q.y = y;

  // normalize bins so that 2D integral gives unity, or
  // a total number of samples in the bin
  q.z = data.bin;
  if (HISTOGRAM_NORM)
  { q.z = q.z ./ dx ./ dy ./ sum(sum(data.bin)); }

  return q;
};

