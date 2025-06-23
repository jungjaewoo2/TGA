//----------------------------------------------------------------------
// ord2
//
// syntax: </a;b;c;d/> = ord2(wn, z, "ss")
//         </den;num/> = ord2(wn, z, "tf")
//
// Generate continuous second order system.
// </A;B;C;D/> = ord2(Wn, Z, "ss") returns the A,B,C,D representation of the
// continuous second order system with natural frequency Wn and 
// damping factor Z.
//
// </DEN;NUM/> = ord2(Wn,Z,"tf") returns the polynomial transfer function of
// the second order system.
//
// See also: rmodel and tf2ss.
//
//----------------------------------------------------------------------
ord2 = function(wn, z, option)
{

  if (!exist(option)) {
     error("Must specify output option: tf or ss.");
  } else {
     if (option != "ss" && option != "tf") {
         error("Output option is either tf or ss.");
     }
  }


  if (option="ss") {
    // Generate state space system
    a = [0, wn;-wn, -2*z*wn];
    b = [0; wn];
    c = [1, 0];
    d = 0;
    return <<a=a;b=b;c=c;d=d>>;
  } else {
    num = 1;
    den = [1, 2*z*wn, wn*wn];
    return <<num=num;den=den>>;
  }
};

