//----------------------------------------------------------------------
//
// lab2str
//
// syntax: labels=lab2str(ulab)
//
// Convert space delimited label format to the more efficient string 
// matrix format.
//
//----------------------------------------------------------------------
lab2str=function(ulab)
{
  len = 12; // Label length parameter

  m = ulab.nr;
  n = ulab.nc;

  if (m==1) {
    // Put labels into row-wise format for MATPRINT
    // First remove extra spaces (delimiters)
    delim = (ulab==" ")
    ndx = find([delim,0]&&[0,delim]);
    if (!isempty(ulab)) { ulab[ndx] = []; }

    // Next, determine position of delimiters
    upos = [0,find(ulab==" "),length(ulab)+1];

    // Finally form row-wise string matrices
    labels = " "*ones(length(upos)-1,len);
    for (i in 1:length(upos)-1) {
      n = upos[i+1]-upos[i]-1;
      if (n<len) {
        labels[i;len-n+1:len] = ulab[upos[i]+[1:n]]; 
      else 
        labels[i;] = ulab[upos[i]+[1:len]]; 
      }
    }

  else 
    labels = [];
    for (i in 1:m) {
      j = n;
      while (ulab[i;j]!="/" && j>1) { j=j-1; }
      if (ulab[i;j] != "/") { lab = ulab[i;j]; else lab = []; }
      while (ulab[i;j]!=" " && j<ns){ j=j+1; lab = [lab,ulab[i;j]]; }
      labels = [labels;[" "*ones(1,len-length(lab)),lab]];
    }

  }
  return labels;
};

