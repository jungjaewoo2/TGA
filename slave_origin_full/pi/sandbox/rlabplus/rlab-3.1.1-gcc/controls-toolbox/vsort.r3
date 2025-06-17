//---------------------------------------------------------------------
//
// vsort
//
// Syntax: </pind;v2/> = vsort(v1,v2,sp)
//
// matches two complex vectors v1 and v2, returning vs2
// which consists of the elements of v2 sorted so that
// they correspond to the elements in V1 in a least squares sense.
// sp is used to test a quick sort method and is equal to
// sp = sum([1:length(indr)].^2); pind=1 is returned if
// a slow sort method has to be applied.
//
//---------------------------------------------------------------------
vsort = function (v1,v2,sp)
{
  pind = 0;
  if (nargs < 3)
  {
     sp = sum([1:length(v1)].^2);
  }
  // Quick Sort
  p = length(v2);
  vones = ones(p,1);
  tmp = abs(vones*v2.'-v1*vones');
  dum = min(tmp);
  indr1 = mini(tmp);
  indr[indr1] = [1:p];

  // Long (accurate) sort
  if (indr*indr' != sp)
  {
     tmp = sort(abs(v2));
     dum = tmp.val;
     jnx = tmp.idx;
     pind = 1;
     for (j in jnx)
     {
         tmp = abs(v2[j]-v1);
         dum = min(tmp);
         inx = mini(tmp);
         indr[inx] = j;
         v1[inx] = 1e30;
     }
  }
  v2 = v2[indr];

  return <<v2=v2; pind=pind>>;
};

