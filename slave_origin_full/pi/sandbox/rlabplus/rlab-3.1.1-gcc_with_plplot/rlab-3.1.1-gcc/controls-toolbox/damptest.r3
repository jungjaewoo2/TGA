//-----------------------------------------------------------------------------
// ****************
// * DAMP Testing *
// ****************
//

rfile feedback damp step

numGs = 1;
denGs = [1,0,0];

numDs = .81*[1,.2];
denDs = [1,2];
tmp   = feedback(conv(numGs,numDs), conv(denGs,denDs),1,1);
numCL = tmp.num;
denCL = tmp.den;
step(numCL,denCL);
damp(denCL);

