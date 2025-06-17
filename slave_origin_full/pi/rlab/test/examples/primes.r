//
// An example that finds all primes less than limit
//

primes = function (limit)
{
  i = 1; j = 0; cnt = 0;
  for(k in 2:limit)
  {
    j = 2;
    while(mod(k,j) != 0)
    {
      j++;
    }
    if(j == k)     // Found prime
    {
      cnt++;
      prime[1;i] = k;
      i++;
    }
  }
  return prime;
};
