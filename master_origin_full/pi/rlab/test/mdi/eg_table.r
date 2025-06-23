//
//
//

NC = 4;
NR = 4;

SR = [2, 2, 2, 2];
SC = [2, 2, 2, 2];

z  = zeros(NR, NC);
s0 = SC;
r0 = shuffle(1:NR, NR);
for (i in r0)
{
  if (sum(s0) == 0)
  { break; }

  s1 = find(s0 > 0);
  "s1 = "
  s1

  c1 = shuffle(1:s1.n, SR[i])';

  "c1 = "
  c1

  z[i; s1[c1]] = ones(z[i; s1[c1]]);
  s0[s1[c1]] = s0[s1[c1]] - 1;

  "s0 = "
  s0

//   pause();
}

"z = "
z

