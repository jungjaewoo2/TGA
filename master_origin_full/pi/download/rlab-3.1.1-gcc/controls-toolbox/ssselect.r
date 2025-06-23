//-----------------------------------------------------------------------------
//
// ssselect
//
// Syntax: A = ssselect(a,b,c,d,inputs,outputs,states)
//         </ae;be;ce;de/> = ssselect(a,b,c,d,inputs,outputs,states)
//
// This routine extracts a subsystem from a larger system. Calling this
// routine as A=ssselect(a,b,c,d,inputs,outputs) will return the state
// space subsystem with the specified inputs and outputs (vectors that
// are input to the routine). Notice that the vectors input and output
// contain indices into the system inputs and outputs.
//
// Calling the routine as A=ssselect(a,b,c,d,inputs,outputs,states) will return
// the state space subsystem with the specified inputs, outputs, and states.
//	and states.
//
// Note: The matrices ae, be, ce, and de are returned in a list.
//
// A.ae = ae matrix.
// A.be = be matrix.
// A.ce = ce matrix.
// A.de = de matrix.
//
// Copyright (C), by Jeffrey B. Layton, 1994
// Version JBL 940405
//-----------------------------------------------------------------------------

require abcdchk 

ssselect = function(a,b,c,d,inputs,outputs,states)
{

   // Check a,b,c,d system to be sure it is valid
   if (nargs < 6 || nargs > 7) {
      error("ssselect: need 6 or 7 arguments");
   }
   msg=abcdchk(a,b,c,d);
   if (msg != "") {
       error("ssselect: "+msg);
   }

   // Create vectors describing new inputs, outputs, and states
   if (nargs == 6) {
       innew=inputs;
       outnew=outputs;
       statenew=1:a.nr;
   }
   if (nargs == 7) {
       innew=inputs;
       outnew=outputs;
       statenew=states;
   }

   // Select system
   if (!isempty(a)) { ae=a[statenew;statenew]; else ae=[]; }
   if (!isempty(b)) { be=b[statenew;innew];    else be=[]; }
   if (!isempty(c)) { ce=c[outnew;statenew];   else ce=[]; }
   if (!isempty(d)) { de=d[outnew;innew];      else de=[]; }

   return << ae=ae; be=be; ce=ce; de=de >>;
};

