//-----------------------------------------------------------------------------
// ********************
// * SSSELECT Testing *
// ********************
//
// Version JBL 950106
//

// =============================================
//                     Test 1:
// =============================================
rfile banner ssselect
printf(" \n");
banner("Test1");
printf(" \n");

A=[33,2,5;23,200,2;9,2,45];
B=[4,5;12,5;82,1];
C=[34,56,2;6,2,112];
D=[2,0;0,19];

inputs=1;
outputs=1;
// Exact solution
aex=A;
bex=[4;12;82];
cex=[34,56,2];
dex=[2];
AT=ssselect(A,B,C,D,inputs,outputs);

// Check the error
printf("selected A matrix check (this matrix should be ZERO)\n");
aex-AT.ae
printf(" \n");

printf("selected B matrix check (this matrix should be ZERO)\n");
bex-AT.be
printf(" \n");

printf("selected C matrix check (this matrix should be ZERO)\n");
cex-AT.ce
printf(" \n");

printf("selected D matrix check (this matrix should be ZERO)\n");
dex-AT.de
printf(" \n");

// =============================================
//                     Test 2:
// =============================================
printf(" \n");
banner("Test2");
printf(" \n");

A=[33,2,5;23,200,2;9,2,45];
B=[4,5;12,5;82,1];
C=[34,56,2;6,2,112];
D=[2,0;0,19];

inputs=1;
outputs=2;
// Exact solution
aex=A;
bex=[4;12;82];
cex=[6,2,112];
dex=[0];
AT=ssselect(A,B,C,D,inputs,outputs);

// Check the error
printf("selected A matrix check (this matrix should be ZERO)\n");
aex-AT.ae
printf(" \n");

printf("selected B matrix check (this matrix should be ZERO)\n");
bex-AT.be
printf(" \n");

printf("selected C matrix check (this matrix should be ZERO)\n");
cex-AT.ce
printf(" \n");

printf("selected D matrix check (this matrix should be ZERO)\n");
dex-AT.de
printf(" \n");

// =============================================
//                     Test 3:
// =============================================
printf(" \n");
banner("Test3");
printf(" \n");

A=[33,2,5;23,200,2;9,2,45];
B=[4,5;12,5;82,1];
C=[34,56,2;6,2,112];
D=[2,0;0,19];

inputs=1;
outputs=2;
states=[2,3];
// Exact solution
aex=[200,2;2,45];
bex=[12;82];
cex=[2,112];
dex=[0];
AT=ssselect(A,B,C,D,inputs,outputs,states);

// Check the error
printf("selected A matrix check (this matrix should be ZERO)\n");
aex-AT.ae
printf(" \n");

printf("selected B matrix check (this matrix should be ZERO)\n");
bex-AT.be
printf(" \n");

printf("selected C matrix check (this matrix should be ZERO)\n");
cex-AT.ce
printf(" \n");

printf("selected D matrix check (this matrix should be ZERO)\n");
dex-AT.de
printf(" \n");


