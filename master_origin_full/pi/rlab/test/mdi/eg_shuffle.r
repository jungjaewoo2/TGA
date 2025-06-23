//
//
//

N = 5;  // total lenght of the set 1:N
k = 3;  // how many indices we need
i = 2;  // exclude this one

// init
sample_norep_1 = function(N,i,k)
{
  r = zeros(1,k+1); // contains entries in order in which they are drawn
  s = zeros(1,k+1); // contains the sorted entries from 'r'
  r[1] = i;
  s[1] = i;

  for (j in 1:k)
  {
    u1  = int((N - j) * uniform()) + 1;

    r[j+1] = u1;

    for (l in [0:j-1])
    {
      r[j+1] = r[j+1] + (r[j+1] >= s[l+1]);
    }

    // array s is sorted for 1:j
    s[j+1] = r[j+1];

    for (l in [j:1:-1])
    {
      if (s[l+1]<s[l])
      {
        // swap the two
        c = s[l+1];
        s[l+1] = s[l];
        s[l] = c;
        else
          break; // no nead to go through the entire 's'
      }
    }
  }

  return r;
};

M = 10;

x = zeros(M,3);
tic();
for (j in 1:M)
{
  x[j;] = sample_norep_1(N,i,k)[2:4];
}
toc()

// now let's sort it:
y = x;
a = max(max(y)) + 1;
a = (a .^ [(y.nc-1): 0 : -1])';
s = y * a;
y = y[sort(s).ind;];



