//
//
//

N = 10;

x = rand(1,N)
sx = sort(x);
sx.val

nx = x;
i3 = shuffle(1:N,3);
nx[i3] = nan(i3);
nx
snx = sort(nx);
snx.val

M=10;
x2 = rand(M,M);
for (i in 1:M)
{
  nn = int(0.5 * M * uniform()) + 1;
  idx_nan = shuffle(1:M,nn)';
  x2[idx_nan;i] = nan(idx_nan);
}
x2
sx2 = sort(x2);
sx2.val



