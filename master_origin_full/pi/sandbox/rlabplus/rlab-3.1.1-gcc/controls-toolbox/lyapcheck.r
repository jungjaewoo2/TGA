
require lyap2 lyap

A=rand(4,4);

E=eig(A);
Eval=[-4.0,-3.0,-2.0,-5.0];
Anew=conj(E.vec')*diag(Eval)*E.vec;
A=real(Anew);
Aold=A;

B=rand(4,4);
E=eig(B);
Eval=[-12,-3.4,-5.8,-1.3];
Bnew=conj(E.vec')*diag(Eval)*E.vec;
B=real(Bnew);
Bold=B;

C=rand(4,4);
E=eig(C);
Eval=[-0.4,-3.4,-9.5,-6.3];
Cnew=conj(E.vec')*diag(Eval)*E.vec;
C=real(Cnew);
Cold=C;

#####################################################################
printf ("\nStart A,B checks\n");
Xab2=lyap2(A,B);
printf ("The following checks ought to be TRUE (1)\n");
Bcheck = all (all (Bold == B))	# these ought to be identical (TRUE)
Acheck = all( all (Aold == A))

printf ("lyap2 solution self-check (this matrix ought to be SMALL\n");
(A*Xab2 + Xab2*A') + B

printf ("The following checks ought to be TRUE (1)\n");
Xab=lyap(A,,B);
Bcheck = all (all (Bold == B))
Acheck = all (all (Aold == A))

printf ("lyap solution self-check (this matrix ought to be SMALL\n");
(A*Xab + Xab*A') + B

printf("More A,B Solution checks.\n");
Xnorm=norm(Xab-Xab2)
printf ("\nCheck lyap against lyap2 (X - X2)\n");
Xab-Xab2

#####################################################################
printf ("\nStart A,B,C checks\n");
Xabc2=lyap2(A,B,C);
printf ("lyap2 solution self-check (this matrix ought to be SMALL\n");
(A*Xabc2 + Xabc2*B) + C

Xabc=lyap(A,B,C);
printf ("\nlyap solution self-check (this matrix ought to be SMALL\n");
(A*Xabc + Xabc*B) + C

printf("\nMore A,B,C Solution checks.\n");
Xnorm=norm(Xabc-Xabc2)
printf ("\nCheck lyap against lyap2 (X - X2)\n");
Xabc-Xabc2

Bcheck = all( all (Bold == B))	# these ought to be identical
Acheck = all (all (Aold == A))
Ccheck = all (all (Cold == C))


