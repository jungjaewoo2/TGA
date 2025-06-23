//--------------------------------------------------------------------------
// 
//  cloop
//
//  syntax: </den;num/> = cloop(num,den)
//          </den;num/> = cloop(num,den,sign)
//          </ac;bc;cc;dc/> = cloop(a,b,c,d)
//          </ac;bc;cc;dc/> = cloop(a,b,c,d,sign)
//          </ac;bc;cc;dc/> = cloop(a,b,c,d,outputs,inputs)
//
//
//	Close unity feedback loops.    -->O-->[System]--+--->
//                                         |             |
//                                         +-------------+
//
//	</Ac;Bc;Cc;Dc/> = CLOOP(A,B,C,D,SIGN) produces a state-space model
//	of the closed-loop system obtained by feeding all the outputs of
//	the system to all the inputs.  If SIGN = 1 then positive feedback
//	is used.  If SIGN = -1 then negative feedback is used.  In all 
//	cases, the resulting system has all the inputs and outputs of the
//	original model.
//
//	</Ac;Bc;Cc;Dc/> = CLOOP(A,B,C,D,OUTPUTS,INPUTS) produces the closed
//	loop system obtained by feeding the specified outputs into the 
//	specified outputs.  The vectors OUTPUTS and INPUTS contain indexes
//	into the outputs and inputs of the system respectively.  Positive
//	feedback is assumed. To close with negative feedback, use negative
//	values in the vector INPUTS.
//
//	</DENc;NUMc/> = CLOOP(NUM,DEN,SIGN) produces the SISO closed loop 
//	system in transfer function form obtained by unity feedback with
//	the sign SIGN.
//
//	See also parallel, series, and feedback.
//
//--------------------------------------------------------------------------
require abcdchk tfchk 

cloop=function(a,b,c,d,e,f)
{

  if(nargs < 2 || nargs > 6) { error("Wrong number of arguments."); } 
  // Determine which syntax is being used.
  if (nargs == 2) {
    // T.F. form w/o sign (Assume negative feedback)
    </den;num/> = tfchk(a,b); 
    sgn = -1;
  }
  if (nargs == 3) {
    // Transfer function form with sign
    </den;num/> = tfchk(a,b); 
    sgn = sign(c);
  }

  // Form closed loop T.F. system 
  if (nargs == 2 || nargs == 3) {
    ac = num; 
    bc = den - sgn*num;
    return <<num=ac;den=bc>>;
    
  } else { if (nargs >= 4)  {
    // State space systems
    msg = abcdchk(a,b,c,d);
    if (msg != "") {error(msg);}
    ny = d.nr; 
    nu = d.nc;
    if (nargs == 4) {
      // System w/o sign (Assume negative feedback)
      outputs = 1:ny; 
      inputs  = 1:nu; 
      sgn = -ones(1,length(inputs));
    }
    if (nargs == 5) {
      % State space form with sign
      outputs = 1:ny; 
      inputs  = 1:nu; 
      sgn = sign(e)*ones(1,length(inputs));
    }
    if (nargs == 6) {
      // State space form w/selection vectors
      outputs = e; 
      inputs = abs(f); 
      sgn = sign(f);
    }

    // Form Closed Loop State-space System
    nin  = length(inputs); 
    nout = length(outputs);
    if (nin != nout) {
      error("The number of feedback inputs and outputs are not equal");
    }
    nx = a.nr; 
    na = a.nc;

    // Form feedback column and row, deal with nonzero D22
    S  = [a,b;c,d];
    Bu = S[;[nx+inputs]]; 
    Cy = S[[nx+outputs];];
    if (!isempty(Cy)) {
      Cy[sgn==-1;] = -Cy[sgn==-1;];  // Get ready for negative feedback
      E = eye(nout,nout) - Cy[;[nx+inputs]];    // E=(I-D22)
      Cy = E\Cy;
    
      Sc = S + Bu*Cy;   //Close Loop

      // Extract closed loop system
      ac = Sc[1:nx; 1:nx];
      bc = Sc[1:nx; nx+1:nx+nu];
      cc = Sc[nx+1:nx+ny; 1:nx];
      dc = Sc[nx+1:nx+ny; nx+1:nx+nu];
    } else {
      ac=a; bc=b; cc=c; dc=d;
    }
    return <<ac=ac;bc=bc;cc=cc;dc=dc>>;
  }}

};
