
require lyap2 lyap  lyapmsf

A=rand(4,4);

E=eig(A);
Eval=[-4.0,-3.0,-2.0,-5.0];
Anew=conj(E.vec')*diag(Eval)*E.vec;
A=real(Anew);

B=rand(4,4);
E=eig(B);
Eval=[-12,-3.4,-5.8,-1.3];
Bnew=conj(E.vec')*diag(Eval)*E.vec;
B=real(Bnew);
B=eye(size(B));

C=rand(4,4);
E=eig(C);
Eval=[-0.4,-3.4,-9.5,-6.3];
Cnew=conj(E.vec')*diag(Eval)*E.vec;
C=real(Cnew);
C=eye(size(C));

//--------------------------------------------------------

printf(" \n");
printf("A,B Solution.\n");

X=lyap(A,,B);

X2=lyap2(A,B);

X3=lyapmsf(A,B);

Xnorm=norm(X-X2)
X2norm=norm(X-X3)
A*X3+X3*A'+B;

//--------------------------------------------------------
printf(" \n");
printf("A,B,C Solution.\n");

X=lyap(A,B,C);

X2=lyap2(A,B,C);

X3=lyapmsf(A,B,C);

Xnorm=norm(X-X2)
X2norm=norm(X-X3)
A*X3+X3*B+C;


