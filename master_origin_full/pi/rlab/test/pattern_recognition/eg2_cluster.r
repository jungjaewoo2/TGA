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
z = 2;

t_0 = m + d * uniform(int(q*n),1);

no = n - int(q*n);
o_0  = m + d + z + (100 - (m+d+z)) * uniform(no,1);

l_0 = [t_0;o_0];
minl_0 = min(l_0);
l_0 = l_0 ./ minl_0;

t  = log10(l_0);


cc2 = cluster.knn(t, 2);
if (length(unique(cc2.feature))!=2)
{ stop("horrible internal error!\n"); }

cc2.val
length(cc2.size)

ls1 = l_0[ find(cc2.feature == 1) ];
 s1 =   t[ find(cc2.feature == 1) ];
lm1 = mean(ls1);
 m1 = mean( s1);
ld1 = max(ls1) - min(ls1);
 d1 = max( s1) - min( s1);

ls2 = l_0[ find(cc2.feature == 2) ];
 s2 = t[ find(cc2.feature == 2) ];
lm2 = mean(ls2);
 m2 = mean( s2);
ld2 = max(ls2) - min(ls2);
 d2 = max( s2) - min( s2);

printf("E(q) = %g\n", max(length(s1), length(s2))/(length(s1) + length(s2)));

gnuwin(1);
lbins=0:100:1;
gnulegend(["Original Set", "Cluster 1", "Cluster 2"]);
gnuplot(<<a=hist(l_0,lbins); b=hist(ls1,lbins); c=hist(ls2,lbins)>>);

abs(m1 - m2) < min(d1,d2)


    