//-----------------------------------------------------------------------------
//
// ssdelete
//
// Syntax: A=ssdelete(a,b,c,d,e,input,outputs)
//         </ar;br;cr;dr/> = ssdelete(a,b,c,d,e,input,outputs)
//
// This routine removes inputs, outputs, and state from a state-space system.
//
// Calling this routine as A=ssdelete(a,b,c,d,input,outputs) will remove the
// specified inputs and outputs from the system. The vectors inputs and
// outputs contain the indices into the system inputs and outputs.
//
// Calling the routine as A=ssdelete(a,b,c,d,inputs,outputs,states) will
// also return a state-space model with the specified inputs, outputs, and
// states removed from the system.
//
//      The results are returned in a list:
//      A.ar = selected A matrix
//      A.br = selected B matrix
//      A.cr = selected C matrix
//      A.dr = selected D matrix
//
//-----------------------------------------------------------------------------

require abcdchk 

ssdelete = function(a,b,c,d,e,f,g)
{

   // Check a,b,c,d system to be sure it is valid
   msg=abcdchk(a,b,c,d);
   if (msg != "") {
       estr="ssdelete: "+msg;
       error(estr);
   }

   // Create vectors describing new inputs, outputs, and states
   if (nargs == 6) {
       inputs=e;
       outputs=f;
       states=[];
   }
   if (nargs == 7) {
       inputs=e;
       outputs=f;
       states=g;
   }

   // Remove the specified inputs,outputs and states


   ar=a; br=b; cr=c; dr=d;
   nx=a.nr; na = a.nc;
   ny=d.nr; nu = d.nc;
   if (length(states)!=nx) {
      if (!isempty(a)) { ar[;states]=[]; ar[states;]=[];  } else { ar=[]; }
      if (!isempty(b)) { br[states;]=[]; br[;inputs]=[];  } else { br=[]; }
      if (!isempty(c)) { cr[;states]=[]; cr[outputs;]=[]; } else { cr=[]; }
   } else {
      ar=[]; br=[]; cr=[];
   }
   if ((length(inputs)!=nu)&&(!isempty(d))) {
     dr[;inputs]=[]; dr[outputs;]=[];
   } else {
     dr=[];
   }

   return <<ar=ar;br=br;cr=cr;dr=dr>>;

};

