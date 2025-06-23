//
// utility functions for im library
//

//
// brightness and contrast formulae from
//    https://en.wikipedia.org/wiki/Image_editing#Contrast_change_and_brightening
//
// if (brightness < 0.0)  value = value * ( 1.0 + brightness);
//                  else value = value + ((1 - value) * brightness);
// value = (value - 0.5) * (tan ((contrast + 1) * PI/4) ) + 0.5;
bright_contr_ch = function(v0, b0, c0, qrange)
{
  global(pi);
  if(!exist(qrange))
  { qrange = 255; }

  b = b0 / qrange;
  c = 2 * c0 / qrange;
  value = v0 / qrange;

  if (b < 0.0)
  {
    value *= ( 1.0 + b );
  }
  else
  {
    value += (1.0 - value) * b;
  }

  value = (value - 0.5) * ( tan((c+1.0)*0.25*pi) ) + 0.5;
  value = ifelse(value>1,1,value);
  value = ifelse(value<0,0,value);
  return int(qrange*value);
};


