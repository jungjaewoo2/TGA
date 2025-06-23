//
//
//
gnuwins (1);

// q = 0.9;
q = 0.5 + 0.5 * uniform();
printf("q = %g\n", q);

n = 100;
m = 1;
d = 2;
z = 1;

t_0 = m + d * uniform(int(q*n),1);

no = n - int(q*n);
o_0  = m + d + z + (100 - (m+d+z)) * uniform(no-1,1);
o_0 = [3.8; o_0];

l_0 = [t_0;o_0];

sl_0 = sort(l_0).val;

// clustering algorhytm
t    = log10(sl_0);
cc2  = cluster.knn(t, 2);

y = [];
for (i in 3:length(sl_0))
{
  del1 = sl_0[i] - 2 * sl_0[i-1] + sl_0[1];
  del2 = sl_0[i] - mean(sl_0[1:(i-1)]) - 3 * var(sl_0[1:(i-1)]) .^ 0.5;
  del3 = mean(sl_0[1:(i)]) - mean(sl_0[1:(i-1)]) - var(sl_0[1:(i-1)]) .^ 0.5;
  y = [y; [i, sl_0[i], del1, del2, del3, cc2.feature[i]]];
}
y

// x_thr = sl_0[i];

// ls1 = len(find(sl_0 <= x_thr));
// printf("E(q) = %g\n", ls1 / len(sl_0));



