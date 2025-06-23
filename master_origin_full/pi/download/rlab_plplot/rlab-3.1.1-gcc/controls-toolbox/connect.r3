//--------------------------------------------------------------------------
//
// connect
//
// Syntax: </aa;bb;cc;dd/> = connect(a,b,c,d,q,iu,iy)
//
// Given the block diagram of a system, CONNECT can be used to
// form an (A,B,C,D) state-space model of the system.
//
// Returns the state-space matrices (Ac,Bc,Cc,Dc) of a system given the 
// block diagonal, unconnected (A,B,C,D) matrices and a matrix Q that 
// specifies the interconnections.  The matrix Q has a row for each 
// block, where the first element of each row is the number of the 
// block.  The subsequent elements of each row specify where the 
// block gets its summing inputs, with negative elements used to 
// indicate minus inputs to the summing junction.  For example, if 
// block 7 gets its inputs from the outputs of blocks 2, 15, and 6, 
// and the block 15 input is negative, the 7'th row of Q would be 
// [7, 2, -15, 6].  The vectors INPUTS and OUTPUTS are used to select 
// the final inputs and outputs for (Ac,Bc,Cc,Dc). 
//   
// See also: blkbuild and cloop.
//
//--------------------------------------------------------------------------

connect = function(a,b,c,d,q,iu,iy)
{
  mq = q.nr; nq = q.nc;
  md = d.nr; nd = d.nc; 

  // Form k from q, the feedback matrix such that u = k*y forms the
  // desired connected system.  k is a matrix of zeros and plus or minus ones.

  k = zeros(nd,md);
  // Go through rows of Q
  for (i in 1:mq)
  {
     // Remove zero elements from each row of Q
     qi = q[i;find(q[i;])];
     // Put the zeros and +-ones in K
     if (qi.nc != 1) {
        k[qi[1];abs(qi[2:n])] = sign(qi[2:n]);
     }
  }

  // Use output feedback to form closed loop system
  //	.
  //	x = Ax + Bu
  //	y = Cx + Du      where  u = k*y + Ur 
  //
  bb = b/(eye(nd,nd) - k*d);
  aa = a + bb*k*c;
  t = eye(md,md) - d*k;
  cc = t\c;
  dd = t\d;

  // Select just the outputs and inputs wanted:
  bb = bb[;iu];
  cc = cc[iy;];
  dd = dd[iy;iu];
  
  return <<aa=aa;bb=bb;cc=cc;dd=dd>>;
};

