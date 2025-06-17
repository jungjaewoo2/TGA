//
//
//
plwins (2);

n = 40;
d = 2;

x = 0.25 * gaussian(n,d);
y = 10 + gaussian(n,d);

t = shuffle([x;y]);

// now pao clustering:

c1 = [];
for (i in 1:15)
{
  c1 = [c1; i, cluster.pao(t, i).val.nr];
}

plwin (2);
plimits(,,0,10);
pltitle ("Clustering - Method Pao");
plformat(["with lines"]);
plxlabel("Pao-clustering Distance");
plylabel("Pao-clustering Number of Centers");
plplot( c1 );

d = 5;
cc1 = cluster.pao(t, d.^2 );
m = cc1.val;

plwin (1);
pltitle ("Clustering - Method Pao");
plimits(-1,15,-1,15);
plxlabel("x");
plxlabel("y");
plegend(["Class1", "Class2", "PAO: cluster centers"]);
plformat(["with points pt 1 ps 4 pc 2", "with points pt 1 ps 4 pc 3", "with points pt 1 ps 4 pc 4"]);
plplot( <<a=x;b=y;c=cc1.val>> );

cc2 = cluster.knn(t, 2);
distance(cc2.val, cc1.val)

cc3 = cluster.iso(t, 2);
distance(cc3.val, cc1.val)









