//
//
//
NITER = 1000;

idxf=0;

N = 50;

for (i in 1:NITER)
{
  spinner();
  c1 = 0.05 * uniform(N-1,1);
  c2 = 1 - 0.05 * uniform(1,1);
  d = shuffle([c1;c2]);

  cc  = cluster.knn(d, 2);
  cc2 = cluster.iso(d, 2);

  mm = mini(cc.val);
  cm = d[ find(cc.feature == mm) ];
  if ( (length(cm)!=(N-1)) && (length(cm)!=1) )
  {
    idxf++;
    //stop("HIE\n");
    printf("*");
    idxf++;
//     stop("Failed!\n");
  }
}

if (idxf>0)
{
  printf("\n Have %g failed feature determinations of %g attempts!\n", idxf, NITER);
  stop();
}


