//----------------------------------------------------------------------
//
// series
//
// Syntax: </a;b;c;d/> = series(a1,b1,c1,d1,a2,b2,c2,d2)
//         </a;b;c;d/> = series(a1,b1,c1,d1,a2,b2,c2,d2,outputs1,inputs2)
//         </den;num/> = series(num1,den1,num2,den2)
//
//	Series connection of two systems.  
//
//		u --->[System1]--->[System2]----> y
//
//	</A;B;C;D/> = series(A1,A2,A3,A4,A2,B2,C2,D2) produces an aggregate
//	state-space system consisting of the series connection of systems
//	1 and 2 that connects all the outputs of system 1 connected to 
//	all the inputs of system 2, u2 = y1.  The resulting system has 
//	the inputs of system 1 and the outputs of system 2.
//
//	</A;B;C;D/> = series(A1,A2,A3,A4,A2,B2,C2,D2,OUTPUTS1,INPUTS2) 
//	connects the two system in series such that the outputs of system
//	1 specified by OUTPUTS1 are connected to the inputs of system 2 
//	specified by INPUTS2.  The vectors OUTPUTS1 and INPUTS2 contain 
//	indexes into the output and inputs of system 1 and system 2 
//	respectively.
// 
//	</DEN;NUM/> = series(NUM1,DEN1,NUM2,DEN2) produces the SISO system
//	in transfer function form obtained by connecting the two SISO 
//	transfer function systems in series.
//
//	See also: append, parallel, feedback and cloop.
//
//----------------------------------------------------------------------
require abcdchk append cloop ssselect tfchk

series = function(a1,b1,c1,d1,a2,b2,c2,d2,e,f)
{
  if (nargs < 4 || nargs > 10) { 
     error("Wrong number of input arguments");
  }

  // --- Determine which syntax is being used ---
  if (nargs == 4)  
  {
    // Form Series connection of T.F. system ---
    A = tfchk(a1,b1); 
    B = tfchk(c1,d1);
    a = conv(A.numc,B.numc);
    b = conv(A.denc,B.denc);
    return <<num=a; den=b>>;
    
  else if (((nargs>=5) && (nargs<=7)) || (nargs==9)) 
  {
    error("Wrong number of input arguments.");

  else if ((nargs==8) || (nargs==10))
  {  
    // State space systems 
    msg = abcdchk(a1,b1,c1,d1);
    if(msg!="") { error("System 1 "+msg); }
    msg = abcdchk(a2,b2,c2,d2);
    if(msg!="") { error("System 2 "+msg); }
 
    // Get number of inputs and outputs
    ny1 = d1.nr; 
    nu1 = d1.nc;
    if (nu1 == 0) { dum = b1.nr; nu1 = b1.nc; ny2 = c1.nr; dum = c1.nc; }
    ny2 = d2.nr;
    nu2 = d2.nc;
    if (nu2 == 0) { dum = b2.nr; nu2 = b2.nc; ny2 = c2.nr; dum = c2.nc; }

    if (nargs == 8) {
      // State space systems w/o selection vectors
      inputs1  = [1:nu1];     
      outputs1 = [1:ny1];
      inputs2  = [1:nu2]+nu1; 
      outputs2 = [1:ny2]+ny1; 
    }
    if (nargs == 10) {
      // State space systems with selection vectors
      inputs1  = [1:nu1];     
      outputs1 = e;
      inputs2  = f+nu1;  
      outputs2 = [1:ny2]+ny1;
    }
  
    // Check sizes
    if (length(outputs1)!=length(inputs2)) {
      error("Series connection sizes don't match.");
    }

    // --- Series Connection ---
    </a;b;c;d/> = append(a1,b1,c1,d1,a2,b2,c2,d2);
    </a;b;c;d/> = cloop(a,b,c,d,outputs1,inputs2);
    return ssselect(a,b,c,d,inputs1,outputs2);
  }}}
};

