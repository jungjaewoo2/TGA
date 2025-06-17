//-----------------------------------------------------------------------
//
// cloop
//
// syntax: </den;num/> = cloop(num,den)
//         </den;num/> = cloop(num,den,sign)
//         </ac;bc;cc;dc/> = cloop(a,b,c,d)
//         </ac;bc;cc;dc/> = cloop(a,b,c,d,sign)
//         </ac;bc;cc;dc/> = cloop(a,b,c,d,outputs,inputs)
//
//
//	Close unity feedback loops.    -->O-->[System]--+--->
//                                         |             |
//                                         +-------------+
//
// This routine closes the loop of either a state-space system or
// a transfer function system for either a continuous or discrete system,
// using unity feedback. Three cases of input arguments are possible:
//
// State Space Case    
// ----------------
//
// Calling the routine with G=cloop(a,b,c,d,sign) will produce a
// state-space model of the closed loop system obtained by feeding
// all of the outputs to all of the inputs. If sign = 1 then positive
// feedback is used. If sign = -1 then negative feedback is used.
// In either case, the resulting system will have all of the same inputs
// and outputs of the original model.
//   Note: A list is returned containing the closed loop model.
//         G.ac = Closed Loop A matrix
//         G.bc = Closed Loop B matrix
//         G.cc = Closed Loop C matrix
//         G.dc = Closed Loop D matrix
//
// State Space Case with Specified Inputs and Outputs    
// --------------------------------------------------
//
// Calling the routine with G=cloop(a,b,c,d,outputs,inputs) will
// produce a state-space model of the closed loop system which
// is obtained by feeding back the specified outputs to the
// specified inputs. The vector "outputs" contains indices of the
// desired outputs and the vector "inputs" contains indices of the
// desired inputs. Positive feedback is assumed in this case.
// If negative feedback is desired, then the vector "inputs" needs
// to have negative values.
//   Note: A list is returned containing the closed loop model.
//         G.ac = Closed Loop A matrix
//         G.bc = Closed Loop B matrix
//         G.cc = Closed Loop C matrix
//         G.dc = Closed Loop D matrix
//
// Transfer Function Case    
// ----------------------
//
// Calling the routine with G=cloop(num,den,sign) produces the SISO
// model of the closed loop system obtained by unity feedback with
// the sign "sign" (positive or negative feedback).
//   Note: A list is returned containing the closed loop model.
//         G.num = Closed loop numerator
//         G.den = Closed Loop denominator
//
// Copyright (C), by Jeffrey B. Layton, 1994-95
// Version JBL 950106
//-----------------------------------------------------------------------

require tfchk abcdchk

cloop = function(a,b,c,d,e,f)
{

// Determine which syntax is being used.
// -------------------------------------

   if (nargs == 2) {
       // transfer function without sign on feedback
       A=tfchk(a,b);
       ac=A.numc;
       bc=A.denc-A.numc;
       return <<num=ac; den=bc>>;
   }
   if (nargs == 3) {
       // transfer function with sign on feedback
       A=tfchk(a,b);
       ac=A.numc;
       bc=A.denc-sign(c)*A.numc;
       return <<num=ac; den=bc>>;
   }

  if (nargs >= 4) {
      // State Space model
      msg=abcdchk(a,b);
      if (msg != "") {
          error(msg);
      }
      // define "input" and "outputs"
      if (nargs == 4) {
          outputs=1:c.nr;
          inputs=1:b.nc;
          sgn2=-1;
      }
      if (nargs == 5) {
          outputs=1:c.nr;
          inputs=[1:b.nc];
          sgn2=sign(e);
      }
      if (nargs == 6) {
          outputs=e;
          inputs=abs(f);
          if (sign(f[1]) > 0 ) {
              sgn2=1;
          else
              sgn2=-1;
          }
      }

      // Form Closed Loop State-space System 
      nout=length(outputs);
      if (length(inputs) != length(outputs)) {
          error("CLOOP: The number of feedback inputs and outputs are not equal");
      }

      dtil=d[outputs;];
      dbar=d[;inputs];
      bbar=b[;inputs];
      ctil=c[outputs;];

      dtilbar=dtil[;inputs];
      if (sgn2 > 0) {
          term1=eye(nout,nout) - dtilbar;
      else
          term1=eye(nout,nout) + dtilbar;
      }

      term2=solve(term1',dbar')';
      term3=solve(term1',bbar')';
      dc=d+sgn2*term2*dtil;
      bc=b+sgn2*term3*dtil;
      cc=c+sgn2*term2*ctil;
      ac=a+sgn2*term3*ctil;
      return << ac=ac; bc=bc; cc=cc; dc=dc >>;
   }

};

