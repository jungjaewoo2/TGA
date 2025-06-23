//----------------------------------------------------------------------
//
// feedback
//
// syntax: </a;b;c;d/> = feedback(a1,b1,c1,d1,a2,b2,c2,d2,e,f)
//         </den;num/> = feedback(num1,den1,num2,den2,sign)
//
// Feedback connection of two systems. 
//
//         u-->O-->[System1]---+--->y
//             |               |
//             +---[System2]<--+
//
// </A;B;C;D/> = feedback(A1,B1,C1,D1,A2,B2,C2,D2,SIGN) produces an 
// aggregate state-space system consisting of the feedback connection
// of the two systems 1 and 2.  Typically, system 1 is a plant and 
// system 2 is a compensator.   If SIGN=1 then positive feedback is
// used. If SIGN=-1 then negative feedback is used.  In all cases, 
// the resulting system has the same inputs and outputs as system 1.
//
// </A;B;C;D/> = feedback(A1,B1,C1,D1,A2,B2,C2,D2,INPUTS1,OUTPUTS1) 
// produces the feedback system formed by feeding all the outputs of
// system2 into the inputs of system 1 specified by INPUTS1 and by 
// feeding the outputs of system 2 specified by OUTPUTS1 into all the
// inputs of system 2.  Positive feedback is assumed.  To connect 
// with negative feedback, use negative values in the vector INPUTS1.
//
// </DEN;NUM/> = feedback(NUM1,DEN1,NUM2,DEN2,SIGN) produces the SISO
// closed loop system in transfer function form obtained by 
// connecting the two SISO transfer function systems in feedback 
// with the sign SIGN.  
//
// See also: cloop, parallel, series, and connect.
//
//----------------------------------------------------------------------
require abcdchk append cloop conv ssselect tfchk

feedback = function(a1,b1,c1,d1,a2,b2,c2,d2,e,f)
{

  if (nargs < 4 || nargs > 10) {
     error("Wrong number of arguments");
  }

  // --- Determine which syntax is being used ---
  if (nargs == 4)  
  { 
    // T.F. w/o sign -- Assume negative feedback
    tmp = tfchk(a1,b1); num1 = tmp.numc; den1 = tmp.denc;
    tmp = tfchk(c1,d1); num2 = tmp.numc; den2 = tmp.denc;
    sgn = -1;
  }
  if (nargs == 5)  
  {
    // Transfer function form with sign
    tmp = tfchk(a1,b1); num1 = tmp.numc; den1 = tmp.denc;
    tmp = tfchk(c1,d1); num2 = tmp.numc; den2 = tmp.denc;
    sgn = sign(a2);
  }

  // --- Form Feedback connection of T.F. system ---
  if ((nargs==4)||(nargs==5))
  {
    a = conv(num1,den2);
    b = conv(den1,den2) - sgn*conv(num1,num2);
    return <<num=a; den=b>>;
  } else { if ((nargs==6) || (nargs==7))
  {
    error("Wrong number of input arguments.");

  } else { if (nargs >= 8)  
  { 
    // State space systems
    msg = abcdchk(a1,b1,c1,d1);
    if (msg!="") { error("System 1 "+msg); }    
    msg = abcdchk(a2,b2,c2,d2);
    if (msg!="") { error("System 2 "+msg); }
    ny1 = d1.nr; nu1 = d1.nc; 
    ny2 = d2.nr; nu2 = d2.nc;
    if (nargs == 8) 
    { // systems w/o sign -- assume negative feedback
      inputs1  = -[1:nu1];     
      outputs1 =  [1:ny1];
      inputs2  =  [1:nu2]+nu1; 
      outputs2 =  [1:ny2]+ny1; 
    }
    if (nargs == 9) 
    { // State space systems with sign
      inputs1  = [1:nu1]*sign(e); 
      outputs1 = [1:ny1];
      inputs2  = [1:nu2]+nu1;     
      outputs2 = [1:ny2]+ny1;
    }
    if (nargs == 10) 
    { // State space systems w/selection vectors
      inputs1  = e;           
      outputs1 = f;
      inputs2  = [1:nu2]+nu1; 
      outputs2 = [1:ny2]+ny1;
    }

    // Check sizes
    if ((length(outputs1)!=length(inputs2)) || (length(outputs2)!=length(inputs1)))
    {
       error("Feedback connection sizes don't match."); 
    }

    // Form Closed-Loop Feedback System
    </a;b;c;d/> = append(a1,b1,c1,d1,a2,b2,c2,d2);
    </a;b;c;d/> = cloop(a,b,c,d,[outputs1,outputs2],[inputs2,inputs1]); // close loops
    </a;b;c;d/> = ssselect(a,b,c,d,[1:nu1],[1:ny1]);
    return <<a=a;b=b;c=c;d=d>>;
  }}}
};


